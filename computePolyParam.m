function [J, w, d] = computePolyParam(dt, route, vff, rhoJ, rhoC)

alpha = dt / mean(route.cellLength);
omega_f = vff * rhoC / (rhoJ - rhoC);

J = [[0, 1-alpha*omega_f, alpha*omega_f];
    [0, 1-alpha*omega_f, 0];
    [0, 1, alpha*omega_f];
    [0, 1-alpha*vff, 0];
    [alpha*vff, 1, alpha*omega_f];
    [alpha*vff, 1, 0];
    [alpha*vff, 1-alpha*vff, 0]];

w = [0;
    alpha*omega_f*rhoC;
    -alpha*omega_f*rhoC;
    alpha*vff*rhoC;
    -alpha*omega_f*rhoJ;
    -alpha*vff*rhoC;
    0];

d = [[(rhoJ - rhoC) / rhoC, 1, -rhoJ];
    [1, 0, -rhoC];
    [0, 1, -rhoC]];