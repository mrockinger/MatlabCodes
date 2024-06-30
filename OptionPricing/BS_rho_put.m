function rhoP=BS_rho_put(S0,T,K,sig,r);
% BS_rho_put computes vega for a put option

d1=( log(S0./K) + (r+0.5*(sig^2)).*T )./ (sig.*sqrt(T));
d2=d1-sig.*sqrt(T);
rhoP=-K.*T.*exp(-r.*T).*cdfn(-d2);
