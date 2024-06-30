function St=simhedg(sig);
% SimHedg.m
% simulates price trajectory at 5 minute level. Horizon over which option is simulated is NbD days.
% once trajectory is constructed extract data for time where one wants to hedge
% creates a module that constructs for a given price series a dynamic hedging strategy.
% apply this model in various situations: mistake on the underlying asset, etc...
clc;

S0=100;
K=100;
r=8; r=log(1+r/100); %continuous time interest rate
sig=0.2;
NbDays=60 %days
mu=0.01; % a drift with 10%
Nb5m=3*5; % alternatievs are multiples of 5. E.g. 5*96 =weekly if =1 take all observations, 
       % if =2, means sampling is every 10'. 
       % If =3 -> 15'
       % if every day use 96.
       % ... every week ->5*96 (assuming 5 working days)
       % each day has 8 hours...
crash=1;    % if >, allow for 10% crash simulation in the middle of sample
       
fprintf(1,'Parameters\n');
fprintf(1,'S0.....................................  %8.4f\n', S0);
fprintf(1,'K......................................  %8.4f\n', K);
fprintf(1,'r the continous compounded rate .......  %8.4f\n', r);
fprintf(1,'Number of days before maturity.........  %8.4f\n',NbDays); 
fprintf(1,'Sampling frequency in 5 minute intervals %8.4f\n ',Nb5m);

% simulate a price path
St=Sim_St(S0,sig,NbDays,mu);

if crash>0;
    % Put a -10% crash in the middle of sample
    T=size(St,1);
    T2=floor(T/2);
    St(T2:T)=St(T2:T)*0.9;
end;

% grab daily stuff and make sure the calibration is OK
s=96:96:(NbDays*96);
DSt=[St(1);St(1+s)];
T=size(DSt,1);
rt=100*log(DSt(2:T)./DSt(1:T-1));

fprintf(1,'Verify that annualized volatility of simulated trajectory is about %8.4f ',sig*100);
sqrt(365)*std(rt)


% extract out of St stuff at the right frequency
% time increment dt given that sampling takes place every Nb5m 
dt=Nb5m*5/(8*60*365); 

% extract out of St corresponding elements
NbS=size(St,1);
s=1:Nb5m:NbS;
St=St(s);

% time till maturity before expiration
tau=0:dt:(size(St,1)*dt); tau=tau';
tau=flipud(tau);
tau(size(tau,1))=0.000001;

[St, Callv, hedgev, Bt, alpha1, Vt, dev]=dynhedg_call(St,K,tau,sig,dt,r);


disp('day,st,call, theta, Bt,alpha,vt,deviation');
ResM=[St Callv hedgev Bt alpha1 Vt dev];
ResM

subplot(2,1,1);
plot(St);
title('Underlying price');

subplot(2,1,2);
s=1:size(Callv,1);
plot(s,Callv,s,Vt);
title('value of call and of replicating strategy');


%==============================================================

function [St,Callv,hedgev,Bt,alphav,Vt,dev]=dynhedg_call(St,K,tau,sigma,dt,r);
% dt indicates time increment, tau the time till option expiration
 
%Create some space to hold variables
T       = size(St,1)
hedgev	= zeros(T,1); %place for hedge ratio
Callv	= zeros(T,1); %vector for call price
Bt		= zeros(T,1); %value of risk-free investment
alphav	= zeros(T,1); %number of units of risk-free investment held
Vt		= zeros(T,1); %value of replicating portfolio
dev		= zeros(T,1); %deviation from call price
s		= 0:(T-1); s=s'; 

Bt      = exp( r*dt*s); %here discount 1 rather than compound

for j=1:T

	Callv(j) = BSCall(St(j),K,sigma,tau(j),r);
    hedgev(j) = bscallhr(K,St(j),sigma,tau(j),r);
    
	if j==1;

		alphav(1)=(Callv(1) - hedgev(1)*St(1) ) / Bt(1); % units in risk-free asset
		Vt(1)=Callv(1);
        
	else

		%dynamic trading strategy 
		Vt(j) = alphav(j-1)*Bt(j)+hedgev(j-1)*St(j); % current value of hedging portfolio
		alphav(j) = ( Vt(j)-hedgev(j)*St(j) )/Bt(j);  % amount of bonds

    end	
end

dev = ( Vt-Callv )./Callv;

%****************************************************************

function St=Sim_St(S0,sig,NbD,mu);
% simulates a geometric brownian motion.
% parameters are calibrated for year. 
% Time increment is in 5 minutes interval.
% NbD is number of days.
% each day is supposed to have 8 trading hours, i.e. 8*12=96 five minute increments

N=96*NbD;       % number of 5 minute increments corresponding to NbD days
D=NbD/(365*N);  %Delta=time between two 5 minutes increments
eps=randn(N,1);
Xt = (mu-0.5.*(sig^2)).*D + sig .* sqrt(D).* eps;
Xt=cumsum(Xt);
St=S0*exp(Xt);
St=[S0;St];
