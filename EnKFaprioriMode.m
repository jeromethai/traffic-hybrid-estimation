function  rhoNext = EnKFaprioriMode(rho, percentStateNoise, J, w, d, rhoJ)

[dim, K] = size(rho); % dim = dimension including ghost cells K = number of ensembles
% taking the mean of rho for the noise
stateNoiseCov = percentStateNoise^2*diag(mean(rho, 2).^2); % assuming uncorrelated
stateNoise = mvnrnd(zeros(dim, 1)', stateNoiseCov, K)'; % state noise ensemble
rhoNext = [rho(1, :); zeros(dim - 2, K); rho(dim, :)];

for l=1:K
    m = s2m(rho2s(rho(:, l), d));
    for i = 2:dim-1
        rhoNext(i, l) = J(m(i-1), :) * [rho(i-1, l); rho(i, l); rho(i+1, l)] + ...
            w(m(i-1));
    end
    rhoNext(:,l) = rhoNext(:,l) + stateNoise(:, l);
    rhoNext(:,l) = max(min(rhoNext(:, l), rhoJ), 0);
end
%{
[m,K] = size(rho);% size of the ensemble
% K is the number of ensembles
% taking the mean of rho for the noise
stateNoiseCov = percentStateNoise^2*diag(mean(rho,2).^2); % assuming uncorrelated
stateNoise = mvnrnd(zeros(m,1)',stateNoiseCov,K)'; % state noise ensemble

rhoNext = zeros(m,K);

for l=1:K
    rhoNext(:,l) = godunovSchemeMode(rho(:,l), J, w, d)+stateNoise(:,l);
    rhoNext(:,l) = max(min(rhoNext(:,l),rhoJ),0);
end
end
%}