function [rho, output] = pCTMEnKFMode(rho_initial, nbEnsembles, steps, dt, route, ...
    vff, rhoJ, rhoC)

output = -1;

% initialize the parameters
dim = length(rho_initial); % dimension (including ghost cells)
percent_dev_meas = 0.05; % percentage of dev of the measurements
percent_dev_state = 0.01; % percentage of dev of the state
[J, w, d] = computePolyParam(dt, route, vff, rhoJ, rhoC);

percent_dev_ens = 0.05; % percentage of deviation of the initial values
ensRho = ...
    mvnrnd(rho_initial',...
    percent_dev_ens^2 * eye(dim),...
    nbEnsembles)'; % initial ensemble
rho = zeros(dim, steps + 1);
rho(:, 1) = mean(ensRho, 2);
measSteps = 0;
tic
for j=1:steps
    fprintf('time step: %i - %i\n', j, steps);
    % A priori step
    ensRho = EnKFaprioriMode(ensRho, percent_dev_state, J, w, d, rhoJ);% a priori
    % A posteriori step
    if(mod(j-1,6)==0)
        measSteps = measSteps+1;
        if(size(route.activeSensors{measSteps},1)~=0)
            Hj = route.observationMatrix(...
                route.activeSensors{measSteps},:);
            Hj = [zeros(size(Hj,1),1), Hj, zeros(size(Hj,1),1)];
            measurements = route.densityMeasured(:,measSteps);
            ensRho = EnKF_aposteriori(ensRho,Hj*[0;measurements;0],Hj,rhoJ,percent_dev_meas);% a posteriori
        end
    end
    rho(:, j + 1) = mean(ensRho, 2);
end
toc
