function F = Pools_solver(t,x,xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,nmax,Phi)
%                       (t,x,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff,1)
% 1->T^sum; 2->H^sum; 3->H^s
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
