function [Te,He,Hse,Neq,nmax,Neqcrit] = equilibrio_nmax(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,Phi,nmax)

Fep = (epsilon)/(nu*Rt-1);
fgl = (lambda_s/Gamma) + 1;
Fxi = (1-xi);


F1 =(nmax*(Rt*Fep + 1))-Phi;
F2 = Gamma*(Rt-1)*fgl/Fxi;
F3 = Rt*lambda_s*(1+Fep);
Hse = F1/(F2-F3);
He = ((lambda_s*Hse+nmax)/(Gamma*(Rt-1)))*((nu+epsilon)*Rt-1)/(nu*Rt-1)-Phi/(Gamma*(Rt-1));
Te = -(lambda_s*Hse+nmax)/(Gamma*(nu*Rt-1));


Neq = nu*Gamma*Rt*Te + lambda_s*Hse + lambda_r*He + nmax;
Neqcrit = Neq;