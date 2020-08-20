function N = newInfections(t,X,Gamma,nu,Rt,lambda_s,lambda_r,nmax,eta)
% obs X = [TSum HSum Hs];
idx = t == floor(t);
T = X(1,idx);
Ht = X(2,idx);
Hs = X(3,idx);
nbt0 = lambda_s*Rt*Hs + lambda_r*Rt*Ht;
ne = zeros(size(T));
for i = 1:length(T)
    if nbt0(i) >=nmax
        ne(i) = eta*nmax;
    else
        ne(i) = eta*nbt0(i);
    end
end
N = nu*Gamma*Rt*T + lambda_s*Hs + lambda_r*Ht + ne;
N = N';