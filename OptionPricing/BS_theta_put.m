function thetaP=BS_theta_put(S0,T,K,sig,r)
% BS_theta_put computes theta for a put option

d1=( log(S0./K) + (r+0.5*(sig^2)).*T )./ (sig.*sqrt(T));
d2=d1-sig.*sqrt(T);
I1=0.5*sig.*S0.*pdfn(d1)./sqrt(T);
I2=-K.*r.*exp(-r.*T).*cdfn(-d2);
thetaP=I1+I2;

