function GB2RND()
% here, load the completed FTSE data matrix 
% the chosen date is March 26, 2004.
% estimates a generalized beta following Liu, Shackleton, Taylor, and Xu.
clc

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

    lb=[ 0.0000001; 0.001; 0.001; 0.001 ]; % 1/a, b, p,q
    ub=[  1; 5000; 10; 10]; %upper bounds for parameters
    A=[-1 0 -1 0; ...
        1 0 0 -1];
    b=[-0.1; -0.1];
%    b0=[ 1/400; 4300; 0.2; 0.2]; 
    b0=[ 1/30; 4000; 5; 2]; 
    Param0=[322.3137 4367.8959   0.1123   0.1160;... 
183.9757 4334.1676   0.1265   0.1054;... 
 37.8915 4276.4255   0.6178   0.4226;... 
 57.4484 4310.5384   0.2844   0.2229;... 
74.1992 4070.5913   0.3161   0.1135   ];
b0=Param0(i,:)'; b0(1)=1/b0(1);
    


%    GB2_Obj1(b0,S0,[KC;KP],CPi,r,T,[C; P],id);
    
    options=optimset('Diagnostics','on','Display','iter','MaxFunEvals',2000);

    [beta,Qmin,exitflag,output,lambda,grad,hessian] =...
      fmincon(@GB2_Obj1,b0,...
      A,b,[],[],lb,ub,[],options,...
      S0,[KC;KP],CPi,r,T,[C; P],id); 
     ParamM=[ParamM; [1/beta(1) beta(2:4)']];
    
    % construction of RND
    a=beta(1); b=beta(2); p=beta(3); q=beta(4);
    f=dens_GB2(z,1/a,b,p,q);
 
    Res=[Res f];  
end
plot(z,Res)
GB2RND=Res;
save GB2RND;

disp('The parameters are ');
Np=rows(ParamM);
for i=1:Np
fprintf('%8.4f %8.4f %8.4f %8.4f \n',ParamM(i,:))
end

%**********************************************************

function y=GB2_Obj1(b0,S0,K,CPi,r,T,CP,id)
%
a = 1/b0(1);
b = b0(2);
p = b0(3);
q = b0(4);
NbStrik=8;

Cemp = CP(1:NbStrik);
Pemp = CP(1+NbStrik:2*NbStrik);

KC = K(1:NbStrik);
KP = K(1+NbStrik:2*NbStrik);

Cth = call_GB2(a,b,p,q,S0,KC,CPi,r,T,CP,id);
Pth = put_GB2(a,b,p,q,S0,KP,CPi,r,T,CP,id);

y1 = Cemp - Cth; %calls
y2 = Pemp - Pth; %puts

e = S0>=KC; % in the money options
y = y1.*e + y2.*(1-e);

dist1=y'*y;

% a
% 1/a
% p
% q
%  p+1/a
%  q-1/a

mart = S0 - b*beta(p+1/a,q-1/a)/beta(p,q) * exp(-r*T);
%mart = log(S0) + r*T - log(b)-log(beta(p+1/a,q-1/a))+log(beta(p,q)) - 1;
%mart=0;

y = dist1 + mart^2;    ;

%*************************************************************

function g=dens_GB2(x,a,b,p,q)
% implements the GB2 of Bookstaber and McDonald
g1 = log(a)-a.*p.*log(b);
g2 = -gammaln(p) - gammaln(q) + gammaln(p+q);
g3 = (a.*p-1).*log(x);
g4 = -(p+q).*log( 1+ (x./b).^a ); 
g  = exp(g1+g2+g3+g4);

function g=cdf_GB2(x,a,b,p,q)
z=a.*log(x./b) - log(1+(x./b).^a);
z=exp(z);
g=betacdf(z,p,q);

%==============================================================

function c=call_GB2(a,b,p,q,S0,K,CPi,r,T,CP,id)
NK=size(K,1);
ov=ones(NK,1);
c=S0 .* ( 1-cdf_GB2( K,a.*ov,b.*ov,(p+1./a).*ov,(q-1./a).*ov ) )-...
     K.*exp(-r*T).*( 1-cdf_GB2(K,a.*ov,b.*ov,p.*ov,q.*ov) );

function p=put_GB2(a,b,p,q,S0,K,CPi,r,T,CP,id)
NK=size(K,1);
ov=ones(NK,1);
c=S0 .* ( 1-cdf_GB2( K,a.*ov,b.*ov,(p+1./a).*ov,(q-1./a).*ov ) )-...
     K.*exp(-r*T).*( 1-cdf_GB2(K,a.*ov,b.*ov,p.*ov,q.*ov) );
p=c - S0.*ov + K.*exp(-r*T);
 
%==============================================================