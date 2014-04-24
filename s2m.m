function m = s2m(s)

dim = length(s) + 1;
m = zeros(dim-2, 1);
for i=1:dim-2    
    m(i) = find([(s(i)==1 && s(i+1)==1);
        (s(i)==1 && s(i+1)==2);
        (s(i)==2 && s(i+1)==1);
        (s(i)==2 && s(i+1)==3);
        (s(i)==3 && s(i+1)==1);
        (s(i)==3 && s(i+1)==2);
        (s(i)==3 && s(i+1)==3)]);
end