% Compute BS via numerical integration of the log-normal density
r=0.08;
r=log(1+r);
T=100/365;
sig=0.25;
S0=100;
K=100;


%support of terminal price of underlying
ST=0.1:0.1:300;
[m,n]=size(ST);

% construction of transition probability
zeta=log(S0)+(r-0.5*(sig^2))*T;
x=(log(ST)-zeta)./(sig*sqrt(T));
pST=1./(sqrt(2*pi).*ST.*sig.*sqrt(T)).*exp(-0.5*(x.^2));
plot(ST,pST);
pause;

% payoff
gST=max([ST-K;zeros(1,n)]);
plot(gST);

C=exp(-r*T)*sum(gST.*pST)*(ST(2)-ST(1));
      
fprintf('Price of call option via numerial integration %9.4f ', C);

