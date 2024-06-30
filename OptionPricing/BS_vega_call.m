function vegaC=BS_vega_call(S0,T,K,sig,r)
% BS_theta_call computes vega for a call option

d1=( log(S0./K) + (r+0.5*(sig^2)).*T )./ (sig.*sqrt(T));
vegaC=S0.*pdfn(d1).*sig.*sqrt(T);

