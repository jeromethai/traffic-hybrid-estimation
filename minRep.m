function [alpha, beta] = minRep(s)

% function that computes the minimal-representation

% arguments:
% s (n X 1): column vector of the mode
% s(i) = 1,2,3 for w,l,d

% returns:
% alpha (n-1 X 1): i -> alpha(i)
% alpha(i) = -1,0,1 for unused, H^c_alpha_i, H_alpha_i
% beta (n X 1): i -> beta(i)
% beta(i) = -1,0,1 for unused, H^c_beta_i, H_beta_i

n = length(s) + 1;

alpha = -1 * ones(n-1, 1);
beta = -1 * ones(n, 1);

switch s(1)
    case 1
        alpha(1) = 1;
        beta(2) = 1;
    case 2
        beta(1) = 1;
        beta(2) = 0;
    case 3
        alpha(1) = 0;
        beta(1) = 0;
end

for i = 2:n-1
    switch s(i-1)
        case 1
            switch s(i)
                case 1
                    beta(i+1) = 1;
                case 2
                    %beta(i,j) = 1; %%
                    beta(i+1) = 0;
            end
        case 2
            switch s(i)
                case 1
                    alpha(i) = 1;
                case 3
                    alpha(i) = 0;
                    %beta(i,j) = 0; %%
            end
        case 3
            switch s(i)
                case 1
                    alpha(i) = 1;
                    beta(i+1) = 1;
                case 2
                    beta(i) = 1;
                    beta(i+1) = 0;
                    beta(i-1) = -1;
                case 3
                    alpha(i) = 0;
                    beta(i) = 0;
                    alpha(i-1) = -1;
            end
    end
end