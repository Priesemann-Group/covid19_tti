function L = maxVp(A)

n = length(A(1,:));
[~,D] = eig(A);
l = NaN(n,1);
for i = 1:n
    l(i) = D(i,i);
end
L = max(real(l));