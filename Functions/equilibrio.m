function [Te,He,Hse,Neq,ne,Neqcrit,fhse] = equilibrio(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,Phi,nmax)

Fep = (epsilon*Rt)/(nu*Rt-1);
Fxi = (1+(Gamma/lambda_s))/(1-xi);
Feta = (Rt-1) * (eta*Rt + Fxi)/(eta*Rt + 1);
M = (Fep + Rt - Feta)^-1;

fhse = M/(lambda_s*(eta*Rt+1));
Hse = Phi*fhse;
He = (Phi/(Gamma*(Rt-1))) * (M*(1+Fep) - 1);
Te = (Phi/(Gamma*(1-nu*Rt))) * M;

He2 = (lambda_s/Gamma)*(Gamma/lambda_s+xi)*Hse/(1-xi);
% He-He2
% He = He2
ne = min(nmax,eta*(lambda_s*Rt*Hse + lambda_r*Rt*He));

Neq = nu*Gamma*Rt*Te + lambda_s*Hse + lambda_r*He + ne;
Neqcrit = nmax*(1+eta*Rt)/(eta*Rt*(1-nu*Rt));

