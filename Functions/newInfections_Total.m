function N = newInfections_Total(t,X,Gamma,nu,epsilon,Rt,Phi)
% obs X = [TSum HSum Hs];
idx = t == floor(t);
T = X(1,idx);
Ht = X(2,idx);
if length(Phi)>1
    N = (nu+epsilon)*Gamma*Rt*T + Gamma*Rt*Ht + Phi(idx);
else
    N = (nu+epsilon)*Gamma*Rt*T + Gamma*Rt*Ht + Phi;
end
N = N';