function [modes, numModes] = minDistModes(rho, d, P, threshold)
% function ratio = minDistModes(rho, d, P, threshold)
% function the computes the reduced set of modes from the state rho
% returns:
% modes: (:,j) -> mode_j
% numModes: number of modes

dim = length(rho);
s = rho2s(rho, d);

[minAlpha, minBeta] = minRep(s); % min H-representation
indAlpha = find(minAlpha >= 0); % indices used in the min representation
indBeta = find(minBeta >= 0);

[~, indAlpha, indBeta] = computeRatio(rho, indAlpha, indBeta, d, P, ...
    threshold);
adj = computeAdj(s, indAlpha, indBeta);
s = [s, adj];

numModes = size(s, 2);
modes = zeros(dim-2, numModes);
for j = 1:numModes
    modes(:, j) = s2m(s(:, j));
end
