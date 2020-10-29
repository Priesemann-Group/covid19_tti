function Phi = PHI(t,tact,sigma,Phi0,Phimax)
if length(tact) == 1
    Phi = Phi0 + (Phimax-Phi0)*pdf('Normal',t,tact,sigma);
else
    nphi = length(tact);
    Phi = Phi0;
    for i = 1:nphi
        Phi = Phi + (Phimax-Phi0)*pdf('Normal',t,tact(i),sigma);
    end
end

