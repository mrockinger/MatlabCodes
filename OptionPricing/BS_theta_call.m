function thetaC=BS_theta_call(S0,T,K,sig,r)
% BS_theta_call computes theta for a call option

d1=( log(S0./K) + (r+0.5*(sig^2)).*T )./ (sig.*sqrt(T));
d2=d1-sig.*sqrt(T);
I1=0.5*sig.*S0.*pdfn(d1)./sqrt(T);
I2=K.*r.*exp(-r.*T).*cdfn(d2);
thetaC=I1+I2;
