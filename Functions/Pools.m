function A = Pools(xi,nu,Rt,Gamma,lambda_s,lambda_r,eta,epsilon)

A = zeros(3);
xim = 1-xi;
A(1,:) = [Gamma*(nu*Rt-1)       (1+eta*Rt)*lambda_r                     (1+eta*Rt)*lambda_s];
A(2,:) = [Gamma*epsilon*Rt      Gamma*(Rt-1)-(1+eta*Rt)*lambda_r     -(1+eta*Rt)*lambda_s];
A(3,:) = [xim*Gamma*epsilon*Rt  xim*Gamma*(Rt-eta*Rt*lambda_r/Gamma)    -(xim*eta*Rt*lambda_s+lambda_s+Gamma)];
