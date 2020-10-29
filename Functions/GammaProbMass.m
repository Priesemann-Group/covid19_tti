function [F,x,CDFg] = GammaProbMass(a,b)
n = 12;
x = 0:1:n; x = x';
F = zeros(n+1,1);
CDFg = gamcdf(x,a,b);
for i = 2:n+1
    F(i) = CDFg(i)-CDFg(i-1);
end

