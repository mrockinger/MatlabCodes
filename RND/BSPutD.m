function bsp = BSPutD(s,k,sigma,tau,r,d);
% bsput Black-Scholes call option formula
% computes the price of a Black-Scholes European put option
% k=strike price
% s=price of underlying asset
% sigma=annual volatility
% tau=time to maturity (as a fraction of a year)
% r=interest rate
% d=continuous dividend payout ratio
stau=sqrt(tau);
d1 = (log(s./k) + ( r - d + 0.5.*(sigma.^2)) .* tau)./(sigma.*stau);
d2 = d1 - sigma.*stau;
bsp= k.*exp(-r.*tau).*cdfn(-d2) - s.*cdfn(-d1);

