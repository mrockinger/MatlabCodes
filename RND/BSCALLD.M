function bsc = BSCallD(s,k,sigma,tau,r,d);
% bscall Black-Scholes call option formula
% computes the price of a Black-Scholes European call option
% k=strike price
% s=price of underlying asset
% sigma=annual volatility
% tau=time to maturity (as a fraction of a year)
% r=interest rate
% d=continuous dividend payout ratio
stau = sqrt(tau);
d1   = (log(s./k) + ( r - d + 0.5*(sigma.^2)) .* tau)./(sigma.*stau);
d2   = d1 - sigma.*stau;
bsc  = s.*cdfn(d1) - k.*exp(-r.*tau).*cdfn(d2);

