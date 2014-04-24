function [rho, P] = combining(modeRho, modeP, modeProb, numModes)

% function that does the miximg step

% arguments:
% modeRho: (:,i) -> rho_i state estimate in mode i
% modeP: (:,:,i) -> P_i covariance in mode i
% transition: (i,j) -> pi_ij transition probabilities
% modeProb: j -> mu_j probability of mode j (column vector)
% dim: dimension of the state
% numModes: number of modes

% returns: 
% mixedRho: (:,j) -> rho_j mixed state in mode i
% mixedP: (:,:,j) -> P_j mixed covariance in mode i

% computes mixedProb: (i,j) -> mu_ij

dim = size(modeRho, 1);
rho = modeRho * modeProb;

P = zeros(dim, dim);
for j=1:numModes
    temp = modeRho(:,j) - rho;
    P = P + (modeP(:,:,j) + temp * temp') * modeProb(j);
end