function bsp = BSPut(s,k,sigma,tau,r);
% bsput Black-Scholes call option formula
% computes the price of a Black-Scholes European put option
% k=strike price
% s=price of underlying asset
% sigma=annual volatility
% tau=time to maturity (as a fraction of a year)
% r=interest rate
echo off;
stau=sqrt(tau);
d1=(log(s./(k.*exp(-r.*tau))))./(sigma.*stau)+0.5*sigma.*stau;
d2=d1-sigma.*stau;
bsp=-s.*cdfn(-d1)+k.*exp(-r.*tau).*cdfn(-d2);
echo on;
