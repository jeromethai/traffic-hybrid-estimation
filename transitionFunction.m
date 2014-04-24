function pi = transitionFunction(mode1, mode2)

smoothing = 5;
pi = 1 / (sum(abs(mode1 - mode2) > 0)+ smoothing);