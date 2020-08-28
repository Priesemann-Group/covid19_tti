function N = newInfections_Total(t,X,Gamma,nu,epsilon,Rt,Phi)
% obs X = [TSum HSum Hs];
T = X(1,:);
Ht = X(2,:);
if length(Phi)>1
    N = (nu+epsilon)*Gamma*Rt*T + Gamma*Rt*Ht + Phi;
else
    N = (nu+epsilon)*Gamma*Rt*T + Gamma*Rt*Ht + Phi;
end
N = N';