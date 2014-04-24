function [rho, output] = pCTMEKFModeLL(rho_initial, steps, dt, route, ...
    vff, rhoJ, rhoC)

output = -1;
load('clusterData/rho7to19K3.mat');
rhoK = rho;

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

output = zeros(steps + 1, 2);

output(1,1) = mvnpdf((rho(2:dim-1,1) - rhoK(2:dim-1,1))', ...
    zeros(1,dim-2), P(2:dim-1,2:dim-1));

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
    % P(2:dim-1,2:dim-1)
    if j == 2863
        P(2:dim-1,2:dim-1) = P(2:dim-1,2:dim-1) + 10^-3 * eye(dim-2);
    end
    output(j+1,1) = mvnpdf((rho(2:dim-1,j+1) - rhoK(2:dim-1,j+1))', ...
        zeros(1,dim-2), P(2:dim-1,2:dim-1));
end
toc

for j=10:(steps+1)
    if output(j,1)>0
        output(1,2) = output(1,2) + log(output(j,1));
    end
end

% K1: 9.7215e+04
% K3: 4.7057e+05
% K5: 7.4672e+05
% K7: 3.9765e+05
% K9: 5.4895e+05
% K11: -4.2554e+05
% K13: 6.9016e+05

%{
K1: [-563700.098313950;]
K2: -inf
K3: [1035484.29513133;]
K5: [997046.023648319;]
K7: [-445664.041497401;]
K11: -820691.081489367
%}

%{
K1: [277076.600335868;]
K2: [343321.427518161;]
K3: [916750.357713981;]
K4: [800914.200597177;]
K5: [901375.049120777;]
K6: [897788.705287965;]
K7: [324291.023062488;]
K8: [868401.880781839;]
K9: [509325.761763882;]
K10: [851334.804759612;]

K1: [2433037.82143942;]
K2: [2838068.54218323;]
K3: [3369086.38055155;]
K4: [3789400.11230024;]
K5: [3888109.82562994;]
K6: [3703403.26789152;]
%}