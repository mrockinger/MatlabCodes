
function VolatSurf()
% in this file load option data and extract the volatility surface 
% the chosen date is March 26, 2004.
clc;

% load data
A=xlsread('FTSE100Mar2604.xls');
%A

%Reshape the data
Ft  = A(2,1); % FT100
Kv  = A(1,2:17); % strikes
Tv  = A(2:6,18)/365; % days till maturity
rtv = log(1+A(2:6,19)/100); % the interest rate in continous time
rtv
%return
CP  = A(2:6,2:17); % all the call and put prices

M   = size(CP,1); % number of dates
N   = size(CP,2); % number of options for one date (both C+P)

Cidx=1:2:N;
Pidx=2:2:N;
Cv=zeros(M*N/2,1);
Pv=zeros(M*N/2,1);
for i=1:5;
    for j=1:N/2;
        Cv(1+(i-1)*N/2+(j-1))=CP(i,1+2*(j-1));
        Pv(1+(i-1)*N/2+(j-1))=CP(i,2+2*(j-1));
    end
end

Opt = [Cv; Pv];

S0  = Ft(1,1)*ones(size(Opt,1),1);
K   = Kv(2:2:N)';

K   = kron(ones(5,1),K);
r   = kron(rtv,ones(size(Cv,1)/5,1));
T   = kron(Tv,ones(size(Cv,1)/5,1));

%AllInfo= [S0 [K; K] Opt [r;r]];
CPind = [ones(size(Cv,1),1) ;zeros(size(Cv,1),1)];
AllInfo= [S0 [K; K] CPind [r; r] [T;T] [Cv; Pv]];

% control of data
for i=1:size(AllInfo,1);
   fprintf('%8.2f %8.2f %3.0f %8.2f %8.2f %8.2f\n', AllInfo(i,:));
end

% some general variables
NbStrik=8;
NbMat  =5;
NbCall =NbMat*NbStrik;
NbPut  =NbMat*NbStrik;
NbOpt  = NbCall+NbPut; % call and puts

% extract for options of different maturities the implicit interest rate and
% dividend ratio by using at the money options

ImplVol=zeros(NbMat,1);
ImplDiv=zeros(NbMat,1);
EstQual=zeros(NbMat,1);
Matv  =zeros(NbMat,1); % for output
for i=1:NbMat;
% find the at the money call and put
% grab first all the call and put options for a given maturity
    AllCall=AllInfo(1+(i-1)*NbStrik:i*NbStrik,:);
    AllPut =AllInfo(NbCall+1+(i-1)*NbStrik:NbCall+i*NbStrik,:);
    Matv(i)=AllCall(1,5);
    
    lb=[ 0.0001; 0.000001 ]; % lower bounds for moments
    ub=[ 1; 1  ]; %upper bounds for moments
    b0=[0.05; 0.1]; %set initial imlied interest to 5% and dividend to 10%
    options=optimset('Diagnostics','on','Display','iter');

    [beta,Qmin,exitflag,output,lambda,grad,hessian] =...
      fmincon(@dist1,b0,...
      [],[],[],[],lb,ub,[],options,...
      AllCall,AllPut); 
    ImplVol(i)=beta(1);        
    ImplDiv(i)=beta(2);
    EstQual(i)=Qmin;
end

% present results
fprintf('-----------------------------------------\n');
fprintf('Days till exp implied volatility and dividends\n');
for i=1:NbMat
    fprintf('%8.4f %8.4f %8.4f\n',365*Matv(i), 100*ImplVol(i),100*ImplDiv(i));
end

subplot(2,1,1);
plot(365*Matv,ImplVol)
title('Term structure of implied volatilities');
subplot(2,1,2);
plot(365*Matv,ImplDiv)
title('Implied dividend for all maturities');

%Shimko's method to compute the implied dividend and interest rate

% for i=1:NbMat;
% % find the at the money call and put
% % grab first all the call and put options for a given maturity
%     C = AllInfo(1+(i-1)*NbStrik:i*NbStrik,6);
%     P = AllInfo(NbCall+1+(i-1)*NbStrik:NbCall+i*NbStrik,6);
%     S0= AllInfo(1,1);
%     y = C-P;
%     T = rows(y);
%     x = [ones(T,1)*S0 K(1:NbStrik)]
%     res = ols(y,x);
%     prt(res) 
% end
%return

% complete info in AllInfo by adding implied dividends
id   = kron(ImplDiv,ones(NbStrik,1));
AllInfo=[AllInfo [id;id]]

% extract for all options the volatility
% uses implied dividends
NbOpt=size(AllInfo,1);
ImpVol=zeros(NbOpt,1);
OptDis=zeros(NbOpt,1);
for i=1:NbOpt
    % invert option price numerically
    % gets parameters
    S0 = AllInfo(i,1);
    K  = AllInfo(i,2);
    CPi= AllInfo(i,3);
    r  = AllInfo(i,4);
    T  = AllInfo(i,5);
    CP = AllInfo(i,6);
    id = AllInfo(i,7);

    lb=[ 0.0001 ]; % lower bounds for moments
    ub=[ 10 ]; %upper bounds for moments
    sig0=0.2; %set initial volatility to 20% 
    options=optimset('Diagnostics','on','Display','iter');

    [beta,Qmin,exitflag,output,lambda,grad,hessian] =...
      fmincon(@dist,sig0,...
      [],[],[],[],lb,ub,[],options,...
      S0,K,CPi,r,T,CP,id); 
    ImpVol(i)=beta;        
    OptDis(i)=Qmin;
end

[ImpVol OptDis]

% present results
%plot(AllInfo(1:NbStrik,2),ImpVol(1:NbStrik),...
%    AllInfo(NbCall+1:NbCall+NbStrik,2),ImpVol(NbCall+1:NbCall+NbStrik))

% select out of the money volatilities
e=AllInfo(:,1)>AllInfo(:,2); % =1 => underlying above K. Call option in the money!

% plot the smiles
Res_C=[]; 
Res_P=[];
Res_im=[];
Res_all=[];
for i=1:NbMat
    idx     = NbStrik*(i-1)+1:NbStrik*(i-1)+NbStrik; idx=idx';
    ivc_im  = ImpVol(idx).*e(idx); %in money option
    ivc_all = ImpVol(idx);
    
    idx     = NbCall+NbStrik*(i-1)+1:NbCall+NbStrik*(i-1)+NbStrik; idx=idx';
    ivp_im  = ImpVol(idx).*(1-e(idx)); % in money put option
    ivp_all = ImpVol(idx);
    
    Res_C=[Res_C ivc_all];
    Res_P=[Res_P ivp_all];
    Res_im=[Res_im (ivc_im + ivp_im)/2];
    Res_all=[Res_all (ivc_all + ivp_all)/2]; % average over all puts and calls
end
size(Res_all)
size(AllInfo(1:NbStrik,2))
size(Matv)

% present implied vols
fprintf('-----------------------------------------\n');
fprintf('Implied volatilities \n');
for i=1:NbMat
    fprintf('%8.4f %8.4f %8.4f %8.4f %8.4f\n', Res_im(i,:));
end
fprintf('-----------------------------------------\n');

cla;
surf(AllInfo(1:NbStrik,2),Matv,Res_C'); rotate3d
surf(AllInfo(1:NbStrik,2),Matv,Res_P'); rotate3d
surf(AllInfo(1:NbStrik,2),Matv,Res_im'); rotate3d
%surf(AllInfo(1:NbStrik,2),Matv,Res_all'); rotate3d

AllInfo1=[AllInfo ImpVol];
save AllInfo1

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function x=dist1(b,Cinfo,Pinfo)
% function dist1 extracts the at the money options and sets up the 
% estimation of the implied volatity and implied dividend
sig=b(1);
d  =b(2);

S0C=Cinfo(:,1);
KC =Cinfo(:,2);
TC =Cinfo(:,5);
rC =Cinfo(:,4);
VC =Cinfo(:,6);

S0P=Pinfo(:,1);
KP =Pinfo(:,2);
TP =Pinfo(:,5);
rP =Pinfo(:,4);
VP =Pinfo(:,6);

VthC = BSCallD(S0C,KC,sig,TC,rC,d);
VthP = BSPutD( S0P,KP,sig,TP,rP,d);

distC = sum((VthC-VC).^2); % euclidian distance between theoretical and actual price
distP = sum((VthP-VP).^2);

x     = distC + distP;

%===============================================================================

function d=dist(sig,S0,K,CPind,r,T,CP,d)
%function dist, defines the distance between the theoretical and actual option price
if CPind==1 %option is a call option
    Vth=BSCallD(S0,K,sig,T,r,d);
else
    Vth=BSPutD(S0,K,sig,T,r,d);
end
d=(Vth-CP)^2; % euclidian distance between theoretical and actual price