% hedge exercice
% computes black-scholes option price
% then implements hedging strategy

K=100
S0=100
r=8
r=log(1+r/100) %continuous time interest rate
T=10 
sigma=0.17


'The price of the call option is : '
bscall(K,S0,sigma,T/365,r)

%Create some space to hold variables
hedgev=zeros(T+1,1); % place for hedge ratio
callv=zeros(T+1,1);  % vector for call price@
Bt=zeros(T+1,1);
St=[100;101;102;101;102;104;102;103;104;105;105];
Day=T+1:-1:0;
Day(T+1)=0.00001;
Day=Day';
alpha=zeros(T+1,1);
Vt=zeros(T+1,1);     %value of replicating portfolio
dev=zeros(T+1,1);    %deviation from call price

j=1;
for j=1:(T+1)
	callv(j)=bscall(St(j),K,sigma,Day(j)/365,r);
	hedgev(j)=bscallhr(K,St(j),sigma,Day(j)/365,r);
	Bt(j)=1*exp(r*(j-1)/365);
	if (j==1);
		alpha(1)=(callv(1)-hedgev(1)*St(1))/Bt(1); %units in risk-free asset
		Vt(1)=callv(1);
	else;                                         %dynamic trading strategy
		Vt(j)=alpha(j-1)*Bt(j)+hedgev(j-1)*St(j);  %current value of hedging portfolio
		alpha(j)=(Vt(j)-hedgev(j)*St(j))/Bt(j);    %amount of bonds
	end;	
end

dev=(Vt-callv)./callv;

j=1:T+1;
plot(j,Vt,'-+',j,callv,'-o');
title('value of call option and of hedging portfolio');

ResM=[callv Vt dev St hedgev alpha Bt];

%fid=fopen('d:/finbox/option/out.out','w'); % for an output file
fid=1;


fprintf(fid,'S0=%8.4f ,K=%8.4f ,r=%8.4f ,sig=%8.4f ,T=%8.4f \n',[S0 K r sigma T]);

fmt1='%12.4f';
fmt =Kron(ones(1,size(ResM,2)),fmt1);
fmt=[fmt '\n'];
fprintf(fid,fmt,ResM');

fid=fclose(fid);