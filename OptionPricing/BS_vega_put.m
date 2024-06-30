function vegaP=BS_vega_put(S0,T,K,sig,r)
% BS_vega_put computes vega for a put option

d1=( log(S0./K) + (r+0.5*(sig^2)).*T )./ (sig.*sqrt(T));
vegaP=S0.*pdfn(d1).*sig.*sqrt(T);