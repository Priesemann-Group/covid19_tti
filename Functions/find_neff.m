function [neff,it] = find_neff(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff0,x0,paso,e,it)

sigma0 = estabilidad_Pools(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff0,x0,1e2*x0(2));
if it == 1
    flag = 1;
    while flag
        neff = neff0-sigma0*paso;                   %neff ser? m?s grande si el sist es inestable
        sigma = estabilidad_Pools(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff,x0,1e2*x0(2));
        if sigma0*sigma <0
            it = it+1;
            flag = 0;
        else
            neff0 = neff;
        end
    end
else
    neff = neff0-sigma0*paso;                   %neff ser? m?s grande si el sist es inestable
    sigma = estabilidad_Pools(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff,x0,1e2*x0(2));
end

neff1 = neff0-sigma0*rand(1)*paso;
sigma1 = estabilidad_Pools(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff1,x0,1e2*x0(2));
if or(sigma1*sigma0<0,sigma*sigma1<0)
    if sigma0
        paso = abs(neff0-neff1);
        neff = neff0;
    elseif sigma1
        paso = abs(neff-neff1);
        neff = neff1;
    end
    if (paso/neff)<=e
        it = it+1;
        return
    else
        disp('loop')
        [neff,it] = find_neff(t,xi,rt,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff1,x0,paso/2,e,it+1);
    end
end