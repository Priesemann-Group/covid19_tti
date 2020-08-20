function x = solver_por_tramos(ti,tf,dt,xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,nmax,Phi,x0)

t = ti:dt:tf;
[~,x] = ode45(@(t,x) Pools_solver(t,x,xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,nmax,Phi),t,x0);
x = x';