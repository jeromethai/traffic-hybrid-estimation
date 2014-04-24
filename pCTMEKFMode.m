function [rho, output] = pCTMEKFMode(rho_initial, steps, dt, route, ...
    vff, rhoJ, rhoC)

output = -1;

% initialize the parameters
dim = length(rho_initial); % dimension (including ghost cells)
percent_dev_state = 0.1; % percentage of std deviation
P = percent_dev_state^2 * eye(dim); P(1, 1)=0; P(dim, dim)=0; % initial covariance matrix
percent_dev_meas = 0.05; % percentage of dev of the measurements
[J, w, d] = computePolyParam(dt, route, vff, rhoJ, rhoC); % poly parameters

rho = zeros(dim, steps + 1);
currentRho = rho_initial;
rho(:, 1) = currentRho;
measSteps = 0;

tic
for j=1:steps
    fprintf('time step: %i - %i\n', j, steps);
    % A priori step
    [currentRho, P] = KFaprioriMode(currentRho, P, ...
        s2m(rho2s(currentRho, d)), percent_dev_state, J, w, rhoJ);% a priori
    % A posteriori step
    if(mod(j-1,6)==0)
        measSteps = measSteps+1;
        if(size(route.activeSensors{measSteps},1)~=0)
            Hj = route.observationMatrix(...
                route.activeSensors{measSteps},:);
            Hj = [zeros(size(Hj,1),1), Hj, zeros(size(Hj,1),1)];
            measurements = route.densityMeasured(:,measSteps);
            [currentRho, P, ~] = KFaposteriori(currentRho, P,...
                Hj*[0;measurements;0] , Hj, rhoJ, percent_dev_meas);% a posteriori
        end
    end
    rho(:, j + 1) = currentRho;
end
toc