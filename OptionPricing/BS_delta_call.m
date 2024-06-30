function deltaC=BS_delta_call(S0,T,K,Sig,r)
% BS_delta_call computes delta for a call option

d1=( log(S0./K) + (r+0.5*(Sig^2)).*T )./ (Sig.*sqrt(T));
deltaC=cdfn(d1);
