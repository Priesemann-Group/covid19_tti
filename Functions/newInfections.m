function N = newInfections(t,X,Gamma,nu,Rt,lambda_s,lambda_r,nmax,eta)
% obs X = [TSum HSum Hs];
T = X(1,:);
Ht = X(2,:);
Hs = X(3,:);
nbt0 = eta*(lambda_s*Rt*Hs + lambda_r*Rt*Ht);
ne = zeros(size(T));
for i = 1:length(nbt0)
    ne(i) = min(nmax,nbt0(i));
end
N = nu*Gamma*Rt*T + lambda_s*Hs + lambda_r*Ht + ne;
N = N';