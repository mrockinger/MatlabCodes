% CallSim1.m
% Parameters are S0, r, sig, T, N
% Valuation of an option using simulation
% allow for antithetic simulations

r=0.08;
r=log(1+r);
T=100/365;
N=1000;
D=T/N;
sig=0.25;
S0=100;
K=100;
MSim=100000; %Number of paths to construct option price

Cv=zeros(MSim,1);

for j=1:MSim

    if (mod(j,1000)==0)
         fprintf('Iteration number %6.0f \n',j);
    end;

    eps=randn(N,1);
    Xt1 = (r-0.5.*(sig.^2)).*D + sig .* sqrt(D).* eps;
    %antithetic random variable
    Xt2 = (r-0.5.*(sig.^2)).*D + sig .* sqrt(D).* (-eps); 
    size(Xt1);

    Xt=[Xt1 Xt2];    
    Nc=size(Xt,2);

    Xt=cumsum(Xt);
    St=S0*exp(Xt);

    payoff=[(St(N,:)-K); zeros(1,Nc)];
    payoff=max(payoff);

    C=mean(payoff)/exp(r*T);


    Cv(j)=C;

end

mCv=mean(Cv);

fprintf('The call option price via simulation %8.4f',mCv);