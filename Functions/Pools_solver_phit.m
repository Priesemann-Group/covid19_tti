function F = Pools_solver_phit(t,x,xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,nmax,Phi0,Phimax,sigma,tact)
%                       (t,x,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff,1)
% 1->T^sum; 2->H^sum; 3->H^s
if nargin<15
    tact = 0;
elseif nargin<14
    sigma = 2;
end
Phi = PHI(t,tact,sigma,Phi0,Phimax);  % default is sigma=2;
F = zeros(3,1);                       
xim = 1-xi;
nbt0 = eta*(lambda_r*Rt*x(2) + lambda_s*Rt*x(3));
if nbt0 >=nmax
    ne = nmax;
else
    ne = nbt0;
end
F(1) = Gamma*(nu*Rt-1)*x(1)         + lambda_r*x(2)                 + lambda_s*x(3)         + ne;
F(2) = Gamma*epsilon*Rt*x(1)        + (Gamma*(Rt-1)-lambda_r)*x(2)	- lambda_s*x(3)         - ne + Phi;
F(3) = xim*Gamma*epsilon*Rt*x(1)    + xim*Gamma*Rt*x(2)             -(lambda_s+lambda_r+Gamma)*x(3)  + xim*(Phi-ne);
