function [rho, output] = pCTMRIMMSimple(rho_initial, steps, dt, route, ...
    vff, rhoJ, rhoC)
% IMM FOR EQUIPROBABLE TRANSITIONS

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
[modesJ, numModes] = findReducedModes(rho_initial, d);
modeRho = repmat(rho_initial, 1, numModes);

modeP = zeros(dim, dim, numModes);
for j = 1:numModes
    modeP(:,:,j) = P;
end

modeProb = ones(numModes, 1) / numModes;

tic
for step = 1:steps
    fprintf('step: %i - %i      modes: %i\n', step, steps, numModes);
    % reduced modes
    fprintf('Reduced modes ...\n');
    modesI = modesJ;
    [modesJ, numModes] = findReducedModes(rho(:, step), d); % adjacent modes
    % mixing step
    fprintf('Mixing ...\n');
    [mixedRho, mixedP] = reducedMixingSimple(modeRho, modeP, modeProb, ...
        modesI, modesJ);
    % A priori step
    fprintf('A priori ...\n');
    modeRho = zeros(dim, numModes);
    modeP = zeros(dim, dim, numModes);
    for j = 1:numModes
        [modeRho(:,j), modeP(:,:,j)] = KFaprioriMode(mixedRho(:,j), ...
            mixedP(:,:,j), modesJ(:,j), percent_dev_state, J, w, rhoJ);% a priori
    end
    % A posteriori step
    fprintf('A posteriori ...\n');
    modeProb = ones(numModes, 1) / numModes;
    if(mod(step-1,6)==0)
        measSteps = measSteps+1;
        if(size(route.activeSensors{measSteps},1)~=0)
            Hj = route.observationMatrix(...
                route.activeSensors{measSteps},:);
            Hj = [zeros(size(Hj,1),1), Hj, zeros(size(Hj,1),1)];
            measurements = route.densityMeasured(:,measSteps);
            modeLikelihood = zeros(numModes, 1);
            for j = 1:numModes
                [modeRho(:,j), modeP(:,:,j), modeLikelihood(j)] = ...
                    KFaposteriori(modeRho(:,j), modeP(:,:,j), ...
                    Hj*[0;measurements;0], Hj, rhoJ, percent_dev_meas);% a posteriori
            end
        end
    modeProb = modeLikelihood / sum(modeLikelihood);  
    end
    % combining step
    fprintf('Combining ...\n');
    [rho(:, step + 1), ~] = combining(modeRho, modeP, modeProb, ...
        numModes);
end
toc