function Shimko_mr()
% implements the RND construction method of Shimko
% set for FTSE data; date is March 26, 2004.
clc;
cla;

% load data
load AllInfo1;

% some general variables
NbStrik=8;
NbMat  =5;
NbCall =NbMat*NbStrik;
NbPut  =NbMat*NbStrik;
NbOpt  = NbCall+NbPut; % call and puts

% extract for options of different maturities the implicit interest rate and
% dividend ratio by using at the money options
z=3800:10:5100; z=z'; %support for RND

e=AllInfo1(:,1)>AllInfo1(:,2);
dist=abs(AllInfo1(:,1)-AllInfo1(:,2));
Res=[];
for i=1:NbMat

    S0 = AllInfo1(1,1);
    KC = AllInfo1(1+NbStrik*(i-1):NbStrik*i,2); %call option
    KP = AllInfo1(NbCall+1+NbStrik*(i-1):NbCall+NbStrik*i,2); %put option
    CPi= AllInfo1(1+NbStrik*(i-1):NbStrik*i,3);
    r  = AllInfo1(1+NbStrik*(i-1),4);
    T  = AllInfo1(1+NbStrik*(i-1),5);
    C = AllInfo1(1+NbStrik*(i-1):NbStrik*i,6); %call option
    P = AllInfo1(NbCall+1+NbStrik*(i-1):NbCall+NbStrik*i,6); %put option
    id = AllInfo1(1+NbStrik*(i-1),7);
    ivc=AllInfo1(1+NbStrik*(i-1):NbStrik*i,8); % implied volatilities
    ivp=AllInfo1(NbCall+1+NbStrik*(i-1):NbCall+NbStrik*i,8); 
    
%    plot([ivc ivp])
     
% get the parameters A0, A1, A2
% get good starting values
e = S0>=KC;
y = ivc.*e + ivp.*(1-e);
Kl= KC.*e + KP.*(1-e);
%plot(Kl,y)

Nobs=rows(y);
x=[ones(Nobs,1) Kl Kl.^2];
res=ols(y,x);
prt(res);
%res.beta
%return

lb = [ -100; -10; 0 ]; % lower bounds for parameters
ub = [  100;  10; 1   ]; %upper bounds for parameters
b0 = res.beta; 
if b0(3)<0
    b0(3)=0.00001;
end
 
 options=optimset('Diagnostics','on','Display','iter','MaxFunEvals',2000,...
      'LargeScale','off','TolX',1e-4,'TolFun',1e-4);
 
 [b,Qmin,exitflag,output,lambda,grad,hessian] =...
              fmincon(@obj_Shimko,b0,...
                      [],[],[],[],lb,ub,[],options,...
                      S0,KC,KP,r,T,C,P,id);
 disp('parameters out of optimization');
 b

 A0=b(1); A1=b(2); A2=b(3);
 sig=A0 + A1*z + A2*z.^2;
 plot(z,sig);
    
  f = shimko_dens2(b,z,S0,KC,KP,r,T,C,P,id); % grab the density
  Res=[Res f]; % store the density
end

plot(z,Res);
ShimkoRND=Res;
save ShimkoRND;

%*******************************************************
    
function y=obj_Shimko(b,S0,KC,KP,r,T,C,P,d);
% compute A0; A1, A2 following Shimko
A0 = b(1);
A1 = b(2);
A2 = b(3);

% generate some option prices below and above given K range
% this is necessary to pin down the poloynomial approximation so that one gets
% positive values
Kg1=(3900:25:4100)'; N1=rows(Kg1);
Kg2=(4850:25:5000)'; N2=rows(Kg2); 

sighC  = A0 + A1*KC + A2*KC.^2;
sighP  = A0 + A1*KP + A2*KP.^2;

sighCe = A0 + A1*Kg1 + A2*Kg1.^2;
sighPe = A0 + A1*Kg2 + A2*Kg2.^2;

e = S0>=KC; %use only in the money options

VthC = BSCallD(S0,KC,sighC,T,r,d); % value of call and put over priced range
VthP = BSPutD( S0,KP,sighP,T,r,d);

VthCe = BSCallD(S0,Kg1,sighCe,T,r,d);% theoretical call and puts
VthPe = BSPutD(S0,Kg2,sighPe,T,r,d);

Ce = BSCallD(S0,Kg1,sighC(1),T,r,d); %entension
Pe = BSPutD( S0,Kg2,sighP(end),T,r,d);

%plot(Kg1,VthCe,Kg1,Ce,KC,VthC,KC,C,KP,P,Kg2,Pe,Kg2,VthPe);
%plot(Kg1,VthCe);title('Kg1,VthCe');pause;
%plot(KC,VthC);title('KC,VthC');pause;
%plot(KC,C);title('KC,C');pause;
%plot(KP,P);title('KP,P');pause;
%plot(Kg1,Ce);title('Kg1,Ce');pause;
%plot(Kg2,Pe);title('Kg2,Pe');pause;
%plot(Kg2,VthPe);title('Kg2,VthPe');pause;

M = S0 - exp(-r*T)*BSCallD(S0,0.001,sighC,T,r,d); %for martingale condition

y = (VthC-C).*e+ (VthP-P).*(1-e);
y = [y; VthCe-Ce;  VthPe-Pe];
y = y'*y + M'*M;

%========================================================

function y=pdfn(x);
% density of a normal with mean 0 and variance 1
y=1/sqrt(2*pi)*exp(-0.5*x.^2);

% ========================================================================


function dens = shimko_dens2(b,z,S0,KC,KP,r,T,C,P,d)
A0  = b(1);  
A1  = b(2);  
A2  = b(3);
sig = A0 + A1*z + A2*z.^2;

st=sqrt(T);
v   = sig*st;;
d1 = (log(S0./z) + (r-d)*T + 0.5*v.^2)./v;
d2 = d1-v;

d1x=-1./(z.*v)+(1-d1./v).*(A1+2*A2*z);
d2x=d1x-(A1+2*A2*z);
corrf=d2x-(A1+2*A2*z).*(1-d2.*d2x)-2*A2*z;
n2=pdfn(d2);
dens=-n2.*corrf;
