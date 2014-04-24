function [ratio, indAlpha, indBeta] = computeRatio(rho, ...
    indAlpha, indBeta, d, P, threshold)

lenAlpha = length(indAlpha);
lenBeta = length(indBeta);

ratio = zeros(lenAlpha + lenBeta, 1);
j = 1;

dRow = d(1, :);
for i = 1:lenAlpha
    ii = indAlpha(i);
    ratio(j) = abs(dRow * [rho(ii); rho(ii+1); 1]) / ...
        (sqrt(dRow(1)^2 + dRow(2)^2) * ...
        (dRow(1)^2 * P(ii, ii) + 2 * dRow(1) * dRow(2) * P(ii, ii+1) + ...
        dRow(2)^2 * P(ii+1, ii+1)));
    j = j + 1;
end

for i = 1:lenBeta
    ii = indBeta(i);
    ratio(j) = abs(rho(ii) + d(2, 3)) / (P(ii, ii));
    j = j + 1;
end

ind = ratio < threshold;
indAlpha = indAlpha(ind(1:lenAlpha));
indBeta = indBeta(ind(lenAlpha+1:end));