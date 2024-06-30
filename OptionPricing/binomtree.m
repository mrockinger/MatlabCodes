function [StM,CtM,CtM11]=BinomTree(S0,K,rho,sig,T,N,CallPut);
% BinomTree.m
% returns the elements corresponding to a binomial tree
% the tree, the call prices along tree, the call price


u=exp(sig*sqrt(T/N));
d=1/u;
Rf=exp(rho*T/N);  % 1+r_f
q=(Rf-d)/(u-d);

% verify parameters
fprintf(1,'verify parameters\n');
fprintf(1,'u= %8.4f d= %8.4f Rf=%8.4f q=%8.4f\n',[u d Rf q]);
fprintf(1,' \n' );

% construction of StM tree
Np1=N+1;
StM=zeros(Np1,Np1);
StM(1,1)=S0;

% fill in diagonal
for i=2:Np1
	StM(i,i)=StM(i-1,i-1)*d;
	i=i+1;
end

% fill in sub-diagonal
for i=1:N
	for j=(i+1):Np1
		StM(i,j)=StM(i,j-1)*u;
    end
end

% compute payoff
if CallPut>0 %have a call option
	payoff=StM(:,Np1)-K;
else %have a put option@
	payoff=K-StM(:,Np1);
end

payoff=[payoff zeros(Np1,1)];
payoff=max(payoff')';

% Initialise tree for call option
CtM=zeros(Np1,Np1);
CtM(:,Np1)=payoff;

% iterate backwards
for j=N:-1:1
	for i=j:-1:1
		CtM(i,j)=(q*CtM(i,j+1)+(1-q)*CtM(i+1,j+1))/Rf;
    end
end

CtM11=CtM(1,1); % the option price