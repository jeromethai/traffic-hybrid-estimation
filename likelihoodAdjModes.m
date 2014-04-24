function [modesJ, numModes] = likelihoodAdjModes(rho, d, modesI, modeProb)
%{
dim = size(modesI, 1) + 2;

%[~, ind] = max(modeProb);
ind = ceil(size(modesI, 2) * rand);
s = m2s(modesI(:, ind));
[minAlpha, minBeta] = minRep(s); % min H-representation
indAlpha = find(minAlpha >= 0); % indices used in the min representation
indBeta = find(minBeta >= 0);
adj = computeAdj(s, indAlpha, indBeta);
s = [s, adj];
s = s(:, randperm(size(s, 2)));
s = s(:, 1:10);
numModes = size(s, 2);
modesJ = zeros(dim-2, numModes);

for j = 1:numModes
    modesJ(:, j) = s2m(s(:, j));
end
%}
