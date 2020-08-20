function [Phi,it] = find_dPhi(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi0,nmax,x0,paso,e,it)

sigma0 = estabilidad_Pools(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi0,nmax,x0,1e2*x0(2));
if it == 1
    flag = 1;
    while flag
        Phi = Phi0+sigma0*paso;                   %neff ser? m?s grande si el sist es inestable
        sigma = estabilidad_Pools(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,nmax,x0,1e2*x0(2));
        if sigma0*sigma <0
            it = it+1;
            flag = 0;
        else
            Phi0 = Phi;
        end
    end
else
    Phi = Phi0+sigma0*paso;                   %neff ser? m?s grande si el sist es inestable
    sigma = estabilidad_Pools(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,nmax,x0,1e2*x0(2));
end

Phi1 = Phi0+sigma0*rand(1)*paso;
sigma1 = estabilidad_Pools(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi1,nmax,x0,1e2*x0(2));
if or(sigma1*sigma0<0,sigma*sigma1<0)
    if sigma0
        paso = abs(Phi0-Phi1);
        Phi = Phi0;
    elseif sigma1
        paso = abs(Phi0-Phi1);
        Phi = Phi1;
    end
    if (paso/Phi)<=e
        it = it+1;
        return
    else
        disp('loop')
        [Phi,it] = find_dPhi(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi1,nmax,x0,paso/2,e,it+1);
    end
end