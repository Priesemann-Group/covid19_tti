function vp = findetacrit(xi,Rt,Gamma,nu,epsilon,lambda_s,x)
A = Pools_Estudio(xi,Rt,Gamma,nu,epsilon,lambda_s,x);
[~,d] = eig(A);
dd = zeros(3,1);
for i = 1:3
    dd(i) = d(i,i);
end
vp = max(dd);