function s = m2s(m)

dim = length(m) + 2;
s = zeros(dim-1, 1);
s(1) = find([(m(1)==1 | m(1)==2);
    (m(1)==3 | m(1)==4);
    (m(1)==5 | m(1)==6 | m(1)==7)]);
for i = 1:dim-2
    s(i+1) = find([(m(i)==1 | m(i)==3 | m(i)==5);
        (m(i)==2 | m(i)==6);
        (m(i)==4 | m(i)==7)]);
end