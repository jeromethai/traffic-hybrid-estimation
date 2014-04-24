function [mixedRho, mixedP] = reducedMixingSimple(modeRho, modeP, modeProb, ...
    modesI, modesJ)

dim = size(modeRho, 1);
numModesI = length(modeProb);
numModesJ = size(modesJ, 2);

% compute mixedRho
mixedRho = sum(repmat(modeProb', dim, 1) .* modeRho, 2) / sum(modeProb);
mixedRho = repmat(mixedRho, 1, numModesJ);

% compute mixedP
P = zeros(dim, dim);
for i = 1:dim
    for j = 1:dim
        for k = 1:numModesI
            P(i, j) = P(i, j) + modeProb(k) * modeP(i, j, k);
        end
    end
end
P = P / sum(modeProb);
mixedP = zeros(dim, dim, numModesJ);
for k = 1:numModesJ
    mixedP(:, :, k) = P;
end