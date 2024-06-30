function bshr = bscallhr(k,s,sigma,tau,r)
% bscallhr hedge ratio of an European call option
echo off;
d1=(log(s./(k.*exp(-r.*tau))))./(sigma.*sqrt(tau))+0.5*sigma.*sqrt(tau);
bshr=cdfn(d1);
echo on;