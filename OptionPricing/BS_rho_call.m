function rhoC=BS_rho_call(S0,T,K,sig,r)
% BS_rho_call computes vega for a call option

d1=( log(S0./K) + (r+0.5*(sig^2)).*T )./ (sig.*sqrt(T));
d2=d1-sig.*sqrt(T);
rhoC=K.*T.*exp(-r.*T).*cdfn(d2);
