function MixRND()
% here, load the completed FTSE data matrix 
% the chosen date is March 26, 2004.
% pursues with a mixture of log-normals
clc;

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

% Now extract densities using mixture of log-normals
Res=[]; %here store the densities
ParamM=[]; % parameters to store
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

    lb=[ -4; -4; 0.0001; 0.0001 ]; % lower bounds for parameters
    ub=[  4; 4; 0.8; 0.8]; %upper bounds for parameters
    b0=[ 0.1; 0.1; 0.4; 0.01]; 
    A=[0 0 -1 1];
    b=0;

    av=0.01:0.1:0.99; av=av'; % start optimizing over a belonging to a grid

    GridRes=[];
    for j=1:size(av,1)
    options=optimset('Diagnostics','on','Display','iter','MaxFunEvals',2000);

    [beta,Qmin,exitflag,output,lambda,grad,hessian] =...
      fmincon(@MD_Obj1,b0,...
      A,b,[],[],lb,ub,[],options,...
      S0,[KC;KP],CPi,r,T,[C; P],id,av(j)); 
  
    GridRes=[GridRes; [av(j) beta' Qmin]];
    end
    
    Nbg=size(GridRes,1);
    for i=1:Nbg
    fprintf('%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f \n',GridRes(i,:));
    end
    
    [mi,miidx]=min(GridRes(:,end));
    b0=GridRes(miidx,1:5); b0=b0';
    
    b0
    
    lb=[ 0.0001; -3; -3; 0.0001; 0.0001 ]; % lower bounds for parameters
    ub=[ 0.9999; 3; 3; 0.9; 0.9];          %upper bounds for parameters
 
    options=optimset('Diagnostics','on','Display','iter','MaxFunEvals',2000);

    [beta,Qmin,exitflag,output,lambda,grad,hessian] =...
      fmincon(@MD_Obj,b0,...
      [],[],[],[],lb,ub,[],options,...
      S0,[KC;KP],CPi,r,T,[C; P],id); 

    ParamM=[ ParamM; beta' ];

    a    = beta(1);
    mu1  = beta(2);
    mu2  = beta(3);
    sig1 = beta(4);
    sig2 = beta(5);
    
    % construction of RND
    f = a * get_LN_RND(z,S0,mu1,T,id,sig1) +...
          (1-a) * get_LN_RND(z,S0,mu2,T,id,sig2);
 
    Res = [ Res f ];  
end
plot(z,Res)
MixRND=Res;
save MixRND;

disp('The parameters are ');
ParamM

%**********************************************************

function y=MD_Obj(b,S0,K,CPi,r,T,CP,id)
%
a      = b(1);
piv    = [a; 1-a];
muv    = b(2:3); muv=muv';
sigv   = b(4:5); sigv=sigv';

NbStrik=8;

Cemp=CP(1:NbStrik);
Pemp=CP(1+NbStrik:2*NbStrik);
KC=K(1:NbStrik);
KP=K(1+NbStrik:2*NbStrik);

Cth=a*BSCallD(S0,KC,sigv(1),T,muv(1),id)+...
        (1-a)*BSCallD(S0,KC,sigv(2),T,muv(2),id);

Pth=a*BSPutD(S0,KP,sigv(1),T,muv(1),id)+...
        (1-a)*BSPutD(S0,KP,sigv(2),T,muv(2),id);
% [CP(1:NbStrik) C1 C2]
% [CP(1+NbStrik:2*NbStrik) P1 P2]
% plot(KC,Cth,'r',KP,Pth,'g',KC,Cemp,'c',KP,Pemp,'y');

y1 = Cemp - Cth; %calls
y2 = Pemp - Pth; %puts

e = S0>=K(1:NbStrik); %as usual in the money options
y = y1.*e + y2.*(1-e);

dist1=y'*y;

mart= S0 - a * BSCallD(S0,0.0001,sigv(1),T,muv(1),id)-...
          (1-a) * BSCallD(S0,0.0001,sigv(2),T,muv(2),id);

y = dist1 + mart^2;    

%===========================================================

function y=MD_Obj1(b,S0,K,CPi,r,T,CP,id,a)
y=MD_Obj([a;b],S0,K,CPi,r,T,CP,id);

%************************************************************
    
function y=get_LN_RND(z,S0,r,T,id,sig);
% the benchmark log-normal density
m = log(S0) + (r-id-0.5*sig^2)*T;
s = sig*sqrt(T);
x = ( log(z)-m )/s;
y=1/s*pdfn(x) ./ z;

function y=pdfn(x);
y=1/sqrt(2*pi)*exp(-0.5*x.^2);

%*************************************************************
