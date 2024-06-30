function SimInsPort(sig);
% SimHedg.m
% simulates price trajectory at 5 minute level. Horizon over which option is simulated is NbD days.
% once trajectory is constructed extract data for time where one wants to hedge
% creates a module that constructs for a given price series a dynamic hedging strategy.
% apply this model in various situations: mistake on the underlying asset, etc...
% Here the data is simulated to display stochastic volatility

S0=3750;
K=S0;
r=3; r=log(1+r/100); %continuous time interest rate
sig=0.25;
mu=0.05; % a drift with 5%
kap=2; % parameters for SV
thet=0.01;
v0=0.01;
rho=-0.5;
psi=0.1;


NbDays=180 %days
Nb5m=96; % alternatievs are multiples of 5. E.g. 5*96 =weekly if =1 take all observations, 
       % if =2, means sampling is every 10'. 
       % If =3 -> 15'... every week ->5*96
crash=-1;    % allow for crash simulation in the middle of sample
fid=fopen('d:/finbox/option/out.out','w');
       
fprintf(fid,'Parameters\n');
fprintf(fid,'S0.....................................  %8.4f\n', S0);
fprintf(fid,'K......................................  %8.4f\n', K);
fprintf(fid,'r the continous compounded rate .......  %8.4f\n', r);
fprintf(fid,'Number of days before maturity.........  %8.4f\n',NbDays); 
fprintf(fid,'Sampling frequency in 5 minute intervals %8.4f\n ',Nb5m);

% simulate a price path
St=Sim_St(S0,sig,NbDays,mu,kap,thet, psi,rho,v0);

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


[St, InsPortv, hedgev, Bt, alpha1, Vt, dev]=dynhedg_InsPort(St,K,tau,sig,dt,r);


NbRep=1:size(St,1); NbRep=NbRep';

ResM=[NbRep tau(1:size(St,1)) St InsPortv hedgev Bt alpha1 Vt dev];

fmt0='%8.4f ';
fmt=Kron(ones(1,size(ResM,2)),fmt0);
fmt=[fmt '\n'];

fprintf(fid,'Nb tau ,st,Insured portfolio, theta, Bt,alpha,vt,deviation \n');
fprintf(fid,fmt,ResM(1:10,:)');

fprintf(fid,'\n\n');
T=size(ResM,1);

fprintf(fid,fmt,ResM(T-9:T,:)');

fid=fclose(fid);

subplot(2,1,1);
plot(St);
title('Underlying price');

subplot(2,1,2);
s=1:size(InsPortv,1);
plot(s,InsPortv,s,Vt);
title('value of insured portfolio and of replicating strategy');


%==============================================================

function [St,InsPortv,hedgev,Bt,alphav,Vt,dev]=dynhedg_InsPort(St,K,tau,sigma,dt,r);
% dt indicates time increment, tau the time till option expiration
 
%Create some space to hold variables
T       = size(St,1)
hedgev	= zeros(T,1); %place for hedge ratio
InsPortv= zeros(T,1); %vector for call price
Bt		= zeros(T,1); %value of risk-free investment
alphav	= zeros(T,1); %number of units of risk-free investment held
Vt		= zeros(T,1); %value of replicating portfolio
dev		= zeros(T,1); %deviation from call price
s		= 0:(T-1); s=s'; 

Bt      = exp( r*dt*s); %here discount 1 rather than compound

for j=1:T

	InsPortv(j) = K*exp(-r*tau(j))+BSCall(St(j),K,sigma,tau(j),r);
    hedgev(j) = bscallhr(K,St(j),sigma,tau(j),r);
    
	if j==1;

		alphav(1)=(InsPortv(1) - hedgev(1)*St(1) ) / Bt(1); % units in risk-free asset
		Vt(1)=InsPortv(1);
        
	else

		%dynamic trading strategy 
		Vt(j) = alphav(j-1)*Bt(j)+hedgev(j-1)*St(j); % current value of hedging portfolio
		alphav(j) = ( Vt(j)-hedgev(j)*St(j) )/Bt(j);  % amount of bonds

    end	
end

dev = 100*( Vt-InsPortv )./InsPortv;

%****************************************************************

function St=Sim_St(S0,sig,NbD ,mu,kap,thet, psi,rho,v0);
% simulates a geometric brownian motion.
% parameters are calibrated for year. 
% Time increment is in 5 minutes interval.
% NbD is number of days.
% each day is supposed to have 8 trading hours, i.e. 8*12=96 five minute increments

N=96*NbD;       % number of 5 minute increments corresponding to NbD days
D=NbD/(365*N);  %Delta=time between two 5 minutes increments
eps1=randn(N,1);
eps2=randn(N,1);
xt=( rho*eps1+sqrt(1-(rho^2))*eps2 )*sqrt(D);
var =zeros(N,1);
var(1,:)=v0;
for j=2:N
    var(j,:)=var(j-1,:) + kap*(thet-var(j-1,:))*D + psi*sqrt(var(j-1,:))*xt(j-1);
end

plot(var);
pause

Xt = (mu-0.5.*var).*D + sqrt(var) .* sqrt(D).* eps1;
Xt=[0;Xt];
Xt=cumsum(Xt);
St=S0*exp(Xt);
