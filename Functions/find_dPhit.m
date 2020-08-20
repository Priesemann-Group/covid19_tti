function [Phimax,it] = find_dPhit(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,nmax,Phimax0,x0,paso,e,it)

sigma0 = estabilidad_Pools_dPhi(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,nmax,x0,Phimax0);
if it == 1
    flag = 1;
    while flag
%         disp('loop1')
        Phimax = Phimax0+sigma0*paso;                   %neff ser? m?s grande si el sist es inestable
        sigma = estabilidad_Pools_dPhi(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,nmax,x0,Phimax);
        if sigma0*sigma <0
            it = it+1;
            flag = 0;
        else
            Phimax0 = Phimax;
        end
    end
else
    Phimax = Phimax0+sigma0*paso;                   %neff ser? m?s grande si el sist es inestable
    sigma = estabilidad_Pools_dPhi(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,nmax,x0,Phimax);
end

Phimax1 = Phimax0+sigma0*rand(1)*paso;
sigma1 = estabilidad_Pools_dPhi(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,nmax,x0,Phimax1);
if or(sigma1*sigma0<0,sigma*sigma1<0)
    if sigma0
        paso = abs(Phimax0-Phimax1);
        Phimax = Phimax0;
    elseif sigma1
        paso = abs(Phi0-Phi1);
        Phimax = Phimax1;
    end
    if (paso/Phi)<=e
        it = it+1;
        return
    else
        disp('loop')
        [Phimax,it] = find_dPhit(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,nmax,Phimax1,x0,paso/2,e,it+1);
    end
end