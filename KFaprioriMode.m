function [rhoNext, Pnext] = KFaprioriMode(rho, P, m, percentStateNoise, ...
    J, w, rhoJ)

dim = length(rho);

% update rho
rhoNext = [rho(1); zeros(dim-2, 1); rho(dim)];
for i=2:dim-1
    rhoNext(i) = J(m(i-1),:) * [rho(i-1); rho(i); rho(i+1)] + w(m(i-1));
end
rhoNext = max(min(rhoNext,rhoJ),0);

% update P
Pnext = zeros(dim, dim);
temp = Pnext;
for i=2:dim-1
    for j=2:dim-1
        temp(i,j) = J(m(i-1),:) * [P(i-1,j); P(i,j); P(i+1,j)];
    end
end
for i=2:dim-1
    for j=2:dim-1
        Pnext(i,j) = [temp(i,j-1), temp(i,j), temp(i,j+1)] * J(m(j-1),:)';
    end
end
Pnext = Pnext + percentStateNoise^2*diag(rho.^2);