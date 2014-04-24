function s = findMode(alpha, beta)

% function that finds the mode from an H-representation
% the H-representation can be in terms of Q_i
% the H-representation can be minimal, or adjacent

% arguments:
% alpha (n-1 X 1): i -> alpha(i)
% alpha(i) = -1,0,1 for unused, H^c_alpha_i, H_alpha_i
% beta (n X 1): i -> beta(i)
% beta(i) = -1,0,1 for unused, H^c_beta_i, H_beta_i

% returns:
% s the mode

n = length(beta);
s = ones(n-1,1);

for j = 1:2
    for i = 1:n-1
        
        a = find([alpha(i)==1 && beta(i)==0;
            beta(i)==1 && beta(i+1)==1;
            alpha(i)==0 && beta(i)==1;
            alpha(i)==1 && beta(i+1)==0;
            beta(i)==0 && beta(i+1)==0;
            alpha(i)==0 && beta(i+1)==1]);
        
        if ~isempty(a)
            switch a
                case 1
                    beta(i+1) = 1;
                case 2
                    alpha(i) = 1;
                case 3
                    beta(i+1) = 0;
                case 4
                    beta(i) = 1;
                case 5
                    alpha(i) = 0;
                case 6
                    beta(i) = 0;
            end
        end
    end
end

for i = 1:n-1
    s(i) = find([alpha(i)==1 && beta(i+1)==1;
        beta(i)==1 && beta(i+1)==0;
        alpha(i)==0 && beta(i)==0]);
end