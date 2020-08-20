function Phi = PHI(t,tact,sigma,Phi0,Phimax)
Phi = Phi0 + (Phimax-Phi0)*pdf('Normal',t,tact,sigma);

