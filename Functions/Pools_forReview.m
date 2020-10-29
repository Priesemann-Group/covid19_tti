function A = Pools_forReview(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon,chi)

A = zeros(4);
xim = 1-xi;
chis = (1-chi*xi)/xim;
A(1,:) = [xim*nu*Rt*Gamma*chis-Gamma     xim*nu*Rt*Gamma*chi         (lambda_r+lambda_s)*(1+xim*eta*Rt*chis)                      xim*Rt*eta*lambda_r*chi];
A(2,:) = [xi*nu*Rt*chis*Gamma            xi*nu*Rt*Gamma*chi-Gamma    (lambda_r+lambda_s)*eta*Rt*xi*chis                           xi*Rt*eta*lambda_r*chi+lambda_r];
A(3,:) = [xim*epsilon*Rt*Gamma*chis      xim*epsilon*Rt*Gamma*chi    -(lambda_r+lambda_s)*(1+xim*eta*Rt*chis)-Gamma+xim*Rt*Gamma*chis  -xim*Rt*eta*lambda_r*chi+xim*Rt*Gamma*chi];
A(4,:) = [xi*epsilon*Rt*Gamma*chis       xi*epsilon*Rt*Gamma*chi     -xi*eta*Rt*(lambda_r+lambda_s)*chis+xi*Rt*Gamma*chis              -lambda_r-xi*Rt*eta*lambda_r*chi+xi*Rt*Gamma*chi-Gamma];