function s = rho2s(rho, d)

dim = length(rho);
s = zeros(dim-1, 1);
for i=1:dim-1
    ind = d * [rho(i); rho(i+1); 1] > 0;
    s(i) = find([ind(1)*ind(3); ind(2)*(1-ind(3)); (1-ind(1))*(1-ind(2))]);
end