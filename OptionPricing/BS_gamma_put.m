function gammaP=BS_gamma_put(S0,T,K,sig,r);
% BS_gamma_put computes gamma for a put option

d1=( log(S0./K) + (r+0.5*(sig^2)).*T )./ (sig.*sqrt(T));
d2=d1-sig.*sqrt(T);
gammaP=pdfn(d1)./(S0.*sig.*sqrt(T));
