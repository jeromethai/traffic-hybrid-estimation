function [rho, output] = pCTMRIMM(rho_initial, steps, dt, route, ...
    vff, rhoJ, rhoC, threshold)


output = -1;

% initialize the parameters
dim = length(rho_initial);
percent_dev_state = 0.1;
P = percent_dev_state^2 * eye(dim); P(1, 1)=0; P(dim, dim)=0;
percent_dev_meas = 0.05; % percentage of dev of the measurements
[J, w, d] = computePolyParam(dt, route, vff, rhoJ, rhoC);

rho = zeros(dim, steps + 1);
rho(:, 1) = rho_initial;
measSteps = 0;

% initial state
numModes = dim;
percent_dev_ens = 0.05; % percentage of deviation of the initial values
modeRho = mvnrnd(rho_initial', percent_dev_ens^2 * eye(dim), numModes)'; % initial ensemble
modesJ = zeros(dim-2, numModes);
for j = 1:numModes
    modesJ(:,j) = s2m(rho2s(modeRho(:,j), d));
end

modeP = zeros(dim, dim, numModes);
for j = 1:numModes
    modeP(:,:,j) = P;
end

modeProb = mvnpdf(modeRho', rho_initial', percent_dev_ens^2 * eye(dim));


output = zeros(steps + 1, 6 + dim - 2);
output(1, 1) = numModes;
output(1, 7:end) = modesJ(:, 1)';

for step = 1:steps
    tic
    fprintf('step: %i - %i      modes: %i\n', step, steps, numModes);
    % reduced modes
    modesI = modesJ;
    [modesJ, numModes] = minDistModes(rho(:, step), d, P, threshold);
    % build transition matrix
    transition = buildTransition(modesI, modesJ);
    % mixing step
    [mixedRho, mixedP] = reducedMixing(modeRho, modeP, modeProb, ...
        modesI, modesJ, transition);
    % A priori step
    modeRho = zeros(dim, numModes);
    modeP = zeros(dim, dim, numModes);
    for j = 1:numModes
        [modeRho(:,j), modeP(:,:,j)] = KFaprioriMode(mixedRho(:,j), ...
            mixedP(:,:,j), modesJ(:,j), percent_dev_state, J, w, rhoJ);% a priori
    end
    % A posteriori step
    modeLikelihood = ones(numModes, 1);
    if(mod(step-1,6)==0)
        measSteps = measSteps+1;
        if(size(route.activeSensors{measSteps},1)~=0)
            Hj = route.observationMatrix(...
                route.activeSensors{measSteps},:);
            Hj = [zeros(size(Hj,1),1), Hj, zeros(size(Hj,1),1)];
            measurements = route.densityMeasured(:,measSteps);
            for j = 1:numModes
                [modeRho(:,j), modeP(:,:,j), modeLikelihood(j)] = ...
                    KFaposteriori(modeRho(:,j), modeP(:,:,j), ...
                    Hj*[0;measurements;0], Hj, rhoJ, percent_dev_meas);% a posteriori
            end
        end
    end
    % update of the mode probability
    modeProb = modeLikelihood .* (transition' * modeProb);
    modeProb = modeProb / sum(modeProb);
    % combining step
    [rho(:, step + 1), P] = combining(modeRho, modeP, modeProb, ...
        numModes);
    output(step + 1, 1) = numModes;
    output(step + 1, 2) = toc;
    output(step + 1, 7:end) = modesJ(:, 1)';
end

ss = rho2s(rho(:, 1), d);
mm = s2m(ss);
aa = rho(:, 1) > rhoC;
output(1, 6) = sum(abs(rho(:, 1) - rhoC) < 0.01);
for step = 1:steps
    s = ss;
    m = mm;
    a = aa;
    ss = rho2s(rho(:, step + 1), d);
    mm = s2m(ss);
    aa = rho(:, step + 1) > rhoC;
    output(step + 1, 3) = sum(abs(ss - s)>0);
    output(step + 1, 4) = sum(abs(mm - m)>0);
    output(step + 1, 5) = sum(abs(aa - a)>0);
    output(step + 1, 6) = sum(abs(rho(:, step + 1) - rhoC) < 0.01);
end