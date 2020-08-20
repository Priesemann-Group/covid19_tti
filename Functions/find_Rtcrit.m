function [rtcrit,it] = find_Rtcrit(t,xi,rt0,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff,x0,paso,e,it)

sigma0 = estabilidad_Pools(t,xi,rt0,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff,x0,1e2*x0(2));
if it == 1
    flag = 1;
    while flag
        rtcrit = rt0+sigma0*paso;                   %neff ser? m?s grande si el sist es inestable
        sigma = estabilidad_Pools(t,xi,rtcrit,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff,x0,1e2*x0(2));
        if sigma0*sigma <0
            it = it+1;
            flag = 0;
        else
            rt0 = rtcrit;
        end
    end
else
    rtcrit = rt0+sigma0*paso;                   %neff ser? m?s grande si el sist es inestable
    sigma = estabilidad_Pools(t,xi,rtcrit,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff,x0,1e2*x0(2));
end

rtcrit1 = rt0+sigma0*rand(1)*paso;
sigma1 = estabilidad_Pools(t,xi,rtcrit1,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff,x0,1e2*x0(2));
if or(sigma1*sigma0<0,sigma*sigma1<0)
    if sigma0
        paso = abs(rt0-rtcrit1);
        rtcrit = rt0;
    elseif sigma1
        paso = abs(rtcrit-rtcrit1);
        rtcrit = rtcrit1;
    end
    if (paso/rtcrit)<=e
        it = it+1;
        return
    else
        disp('loop')
        [rtcrit,it] = find_Rtcrit(t,xi,rtcrit1,Gamma,nu,epsilon,lambda_s,lambda_r,eta,Phi,neff,x0,paso/2,e,it+1);
    end
end