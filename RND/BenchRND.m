function BenchRND()
% here, load the completed FTSE data matrix 
% the chosen date is March 26, 2004.
% computes the simplest log-normal density using at the money options
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

e=AllInfo1(:,1)>AllInfo1(:,2);
dist=abs(AllInfo1(:,1)-AllInfo1(:,2));
Res=[];
for i=1:NbMat
    iv=AllInfo1(:,8);

    idx=NbStrik*(i-1)+1:NbStrik*(i-1)+NbStrik; idx=idx';
    ivc=iv(idx).*e(idx);
    
    idx=NbCall+NbStrik*(i-1)+1:NbCall+NbStrik*(i-1)+NbStrik; idx=idx';
    ivp=iv(idx).*(1-e(idx));
    iv=(ivc+ivp)/2; 
    
    sig=min(iv);
    S0 = AllInfo1(1,1);
    r  = AllInfo1(NbStrik*(i-1)+1,4);
    T  = AllInfo1(NbStrik*(i-1)+1,5);
    id = AllInfo1(NbStrik*(i-1)+1,7);
    f=get_LN_RND(z,S0,r,T,id,sig);
    Res=[Res f];   
end
plot(z,Res);
BenchRND=Res;
save BenchRND;

%*******************************************************
    
function y=get_LN_RND(z,S0,r,T,id,sig);
% the benchmark log-normal density
m = log(S0) + (r-id-0.5*sig^2)*T;
s = sig*sqrt(T);
x = ( log(z)-m )/s;
y=1/s*pdfn(x) ./ z;

function y=pdfn(x);
y=1/sqrt(2*pi)*exp(-0.5*x.^2);

%**********************************************************
