function transition = buildTransition(modesI, modesJ)

numModesI = size(modesI, 2);
numModesJ = size(modesJ, 2);
transition = zeros(numModesI, numModesJ);

for i = 1:numModesI
   for j = 1:numModesJ
       transition(i, j) = transitionFunction(modesI(:, i), modesJ(:, j));
   end
end