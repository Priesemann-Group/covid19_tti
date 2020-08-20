function sigma = estabilidad_Pools_dPhi(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi0,neff,x0,Phimax)
[~,X] = ode45(@(t,x) Pools_solver_phit(t,x,xi,nu,rt,Gamma,lambda_s,lambda_r,eta,epsilon,neff,Phi0,Phimax),t,x0);
m = X(end,1)-X(end-1,1);
if m>1
    sigma = -1;
else
    sigma = 1;
end