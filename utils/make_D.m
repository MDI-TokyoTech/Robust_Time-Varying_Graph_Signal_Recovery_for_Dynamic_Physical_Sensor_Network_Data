function D = make_D(M)

D = zeros(M,M-1);
for i = 1:M-1
    D(i,i) = -1;
    D(i+1,i) = 1;
end
