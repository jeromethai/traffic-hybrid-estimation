function [mixedRho, mixedP] = reducedMixing(modeRho, modeP, modeProb, ...
    modesI, modesJ, transition)

dim = size(modeRho, 1);
numModesI = length(modeProb);
numModesJ = size(modesJ, 2);

% compute mixing probability
mixingProb = transition .* repmat(modeProb, 1, numModesJ);
for j = 1:numModesJ
    mixingProb(:, j) = mixingProb(:, j) / sum(mixingProb(:, j));
end

% compute mixedRho
mixedRho = modeRho * mixingProb;

% compute mixedP
mixedP = zeros(dim, dim, numModesJ);
for j = 1:numModesJ
    P = zeros(dim, dim);
    for i = 1:numModesI
        temp = modeRho(:, i) - mixedRho(:,j);
        P = P + (modeP(:, :, i) + temp * temp') * mixingProb(i, j);
    end
    mixedP(:, :, j) = P;
end