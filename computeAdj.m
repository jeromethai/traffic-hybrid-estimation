function adj = computeAdj(s, indAlpha, indBeta)

[alpha, beta] = findRep(s); % H-representation
lenAlpha = length(indAlpha);
lenBeta = length(indBeta);
numAdj = lenAlpha + lenBeta; % number of adjacent modes
adj = zeros(length(s), numAdj);
j = 1;
% adjacent modes by taking the complementary of half-spaces alpha
for i = 1:lenAlpha
    adjAlpha = alpha;
    adjAlpha(indAlpha(i)) = 1 - alpha(indAlpha(i));
    adj(:,j) = findMode(adjAlpha, beta);
    j = j + 1;
end

% adjacent modes by taking the complementary of half-spaces beta
for i = 1:lenBeta
    adjBeta = beta;
    adjBeta(indBeta(i)) = 1 - beta(indBeta(i));
    adj(:,j) = findMode(alpha, adjBeta);
    j = j + 1;
end