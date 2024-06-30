function mCv=CallSim2(S0,K,sig,r,T,N,MSim);
% computes the price of a call option via simulation.
% code is presented as a function.

D=T/N;
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