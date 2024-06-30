function SemPDen()
% SemPDen.m
%
% computes the risk neutral densities using the approach of Abken, Madan,
% and Ramamurtie (1996) Estimation of rsik-neutral and statistical densities
% by Hermite polynomial approximation: With an application to Eurodollar
% Futures Options, Atlanta FED, mimeo.
%
% Imposes positivity restrictions for Hermite polynomial.
% Allows for a free estimation of mean parameter
% Estimation for 250497 using a free mean m. Version is thus not completely
% This program can be freely distributed, modified, improved, whatever....
% Please, report bugs or improvements to Michael Rockinger (Rockinger@hec.fr)
% If in academic work you mention our paper
% "Estimating Gram-Charlier expansions under positivity constraints"
% Jondeau, E., and M. Rockinger, 1998, HEC working paper
% then it would be nice. The Hermite approximation has been developed by Abken
% Madan, and Milne.
% To run program, must have sk.fmt containing the skewness-kurtosis boundary
%
clc; 

%**********************************************************************
%*                   load data and set global parameters              *
%**********************************************************************
%globals for skewness+kurtosis boundaries
eps=0.0001; % distance from boundary
s6=sqrt(6);
s24=sqrt(24);

load sk;

% check that data is correct
%plot(sk(:,1),sk(:,2));
%sk(1:2,:)
%sk(end-1:end,:)


%computes linear interpolation to improve approximation of skewness bound%
[ai,bi]=init_ab(sk,eps);
ai(1:2)
bi(1:2)
ai(end-1:end)
bi(end-1:end)

%boundaries for kurtosis%
kumin=0; kumax=4/s24;

%*****************************************************************************
% load data
load AllInfo1;

% some general variables
NbStrik= 8;
NbMat  = 5;
NbCall = NbMat*NbStrik;
NbPut  = NbMat*NbStrik;
NbOpt  = NbCall+NbPut; % call and puts

% extract for options of different maturities the implicit interest rate and
% dividend ratio by using at the money options
z=3800:10:5100; z=z'; %support for RND

sempm = []; % here store final risk-neutral densities
Param = []; % here store parameters
for i = 1:1%NbMat;
    
S0 = AllInfo1(1,1);
KC = AllInfo1(1+NbStrik*(i-1):NbStrik*i,2); %call option
KP = AllInfo1(NbCall+1+NbStrik*(i-1):NbCall+NbStrik*i,2); %put option
CPi= AllInfo1(1+NbStrik*(i-1):NbStrik*i,3);
r  = AllInfo1(1+NbStrik*(i-1),4);
T  = AllInfo1(1+NbStrik*(i-1),5);
C  = AllInfo1(1+NbStrik*(i-1):NbStrik*i,6); %call option
P  = AllInfo1(NbCall+1+NbStrik*(i-1):NbCall+NbStrik*i,6); %put option
id = AllInfo1(1+NbStrik*(i-1),7);
% To make further computations easier
% convert put into call and consolidate parameters

C2 = P + S0 - exp(-r*T)*KC;
e  = S0 >= KC; % in the money options
C  = C.*e + C2.*(1-e);
K  = KC;
% plot(K,C);
% return

%To make further computations easier
C = C;%*exp(r*T);

A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
options = optimset('Diagnostics','on','Display','iter');
lb = [-0.1; 0.025; -10; -10 ]; 
ub = [ 0.1; 0.5;    30;   30];
param0=[0.002;   0.170;   -0.559;   -0.695 ];

How_is_fit(param0,kumin,kumax,sk,ai,bi,T,S0,K,C,r)
%pause

[b,fmin,retcode,output,lambda,gr,hessian] =...
      fmincon(@obj,param0,A,b,Aeq,beq,lb,ub,nonlcon,options,...
                kumin,kumax,sk,ai,bi,T,S0,K,C,r);
     Param=[Param; b'];
     
fprintf('Optimal parameter');
b

% verify fit
How_is_fit(b,kumin,kumax,sk,ai,bi,T,S0,K,C,r)

b1=b;
b1(4)=rect_x(b1(4),kumin,kumax);       % map kurtosis into allowed domain%
sku=interk(b1(4),sk,ai,bi);            % get upper bound for skewness%
b1(3)=rect_x(b1(3),-sku,sku);          % restrict skewness%
b1';

m  = b(1) ;
s  = b(2);
b0 = 1;
b1 = 0;
b2 = 0;
b3 = b(3);
b4 = b(4);


% recover the traditional risk neutral density for S_T
mu   = log(S0) + ( m - 0.5 * s^2 ) * T;
si   = s*sqrt(T);
tmp  = ( log(z) - mu )/si;
semp = 1/si * ri_ne_de(tmp,[b0;b1;b2;b3;b4],kumin,kumax,sk,ai,bi);
semp = semp./z;

sempm = [sempm semp];

end

SempRND=sempm;
save SempRND;

plot(z,sempm);
title('Risk Neutral Densities for FTSE data');
xlabel('Level of FTSE');
ylabel('Risk Neutral Density');

r=rows(Param);
for i=1:r
    fprintf('%8.4f %8.4f %8.4f %8.4f \n',Param(i,:));
end

%******************************************************************************

function [ai,bi]=init_ab(sk,eps); %initializes parameters for linear interpolation
N  = rows(sk);
k  = sk(:,1);
s  = sk(:,2) - eps;

dk = k(2:N) - k(1:N-1);
ai = ( s(1:N-1) .* k(2:N) - k(1:N-1) .* s(2:N) ) ./ dk;
bi = ( s(2:N) - s(1:N-1) ) ./ dk;

%-----------------------------------------------------------------------------

function f = obj(b,kumin,kumax,sk,ai,bi,T,S0,K,C,r);

m  = b(1); 
s  = b(2);
p0 = 1; 
p3 = b(3); 
p4 = b(4);

p4  = rect_x( b(4), kumin, kumax);       %map kurtosis into allowed domain
sku = interk(p4,sk,ai,bi);               %get upper bound for skewness
p3  = rect_x( b(3), -sku, sku);          %restrict skewness

[a0,a1,a2,a3,a4] = ak(m,s,T,S0,[0.001;K],r);       %for all exercices prices simultaneously

Cth= p0 * a0 + p3 * a3 + p4 * a4;

dist = [S0;C] - Cth;

f    = dist' * dist;

%-----------------------------------------------------------------------------

function [a0,a1,a2,a3,a4]=ak(m,s,T,S0,K,r);
%returns a0, a1, a2, a3, a4

st = s*sqrt(T);
sp = 1/sqrt(2*pi);
d1 = ( log(S0./K) + ( m + 0.5*(s.^2)).*T )./st;
d2 = d1 - st; 

f   = S0; 
d1f = st.*f; 
d2f = st.*d1f; 
d3f = st.*d2f; 
d4f = st.*d3f;

n1   =  sp.*exp(-0.5*d1.^2); 
d1n1 = -d1.*n1; 
d2n1 = (d1.^2-1).*n1; 
d3n1 = (3*d1 - d1.^3).*n1;

n2   = sp.*exp(-0.5*d2.^2); 
d1n2 = -d2.*n2; 
d2n2 = (d2.^2-1).*n2; 
d3n2 = (3*d2-d2.^3).*n2;

K  = K.*exp(-r.*T);
a0 = f*cdfn(d1) - K.*cdfn(d2);

a1 = d1f.*cdfn(d1) + f.*n1 - K.*n2;

a2 = d2f.*cdfn(d1) + 2*d1f.*n1 + f.*d1n1 - K.*d1n2; 
a2 = a2/sqrt(2);

a3 = d3f.*cdfn(d1) + 3*d2f.*n1 + 3*d1f.*d1n1 + f.*d2n1 - K.*d2n2; 
a3 = a3/sqrt(6);

a4 = d4f.*cdfn(d1) + 4*d3f.*n1 + 6*d2f.*d1n1 + 4*d1f.*d2n1 + f.*d3n1 - K.*d3n2; 
a4 = a4/sqrt(24);

%-----------------------------------------------------------------------------

function Qz=ri_ne_de(z,b,kumin,kumax,sk,ai,bi); % evaluates the density Q
b0 = b(1); 
b1 = b(2); 
b2 = b(3); 
b3 = b(4); 
b4 = b(5);

b4  = rect_x(b(4),kumin,kumax);       %map kurtosis into allowed domain%
sku = interk(b4,sk,ai,bi);                    %get upper bound for skewness%
b3  = rect_x(b(3),-sku,sku);          %restrict skewness%

Qz  = 1/sqrt(2*pi)*exp(-0.5*z.^2).*( b0-b2/sqrt(2)+3*b4/sqrt(24)+...
      (b1-3*b3/sqrt(6)).*z + (b2/sqrt(2)-6*b4/sqrt(24)) * z.^2+...
      b3/sqrt(6).*z.^3 +...
      b4/sqrt(24).*z.^4 );

%-----------------------------------------------------------------------------

function x=interk(k,sk,ai,bi);
%returns skewness (>0) for given k, in (0,4), using linear interpolation%
s6=sqrt(6);
s24=sqrt(24);

if k<0
    disp('procedure does not allow negative kurtosis'); 
    return;
end
if k >= (4/s24)
    x=0.0;
else
[ma,i]=max( k < sk(:,1) );
i=i-1;
x=ai(i)+bi(i)*k;
end

%-----------------------------------------------------------------------------

function y = rect_x(x,a,b)
% x is a vector with elements ranging between real -infty, +infty.
% Perform logistic and map into ]a,b[ 
ex = exp(x);
y  = a + (ex./(1+ex)) * (b-a);

%-----------------------------------------------------------------------------

function How_is_fit(b,kumin,kumax,sk,ai,bi,T,S0,K,C,r)
% traces theoretical and empirical option prices
hold off;
cla;

 S0
 K
 C
 b
 kumin
 kumax
 T

m  = b(1); 
s  = b(2);
p0 = 1; 
p3 = b(3); 
p4 = b(4);

p4  = rect_x( b(4), kumin, kumax);       %map kurtosis into allowed domain
sku = interk(p4,sk,ai,bi);               %get upper bound for skewness
p3  = rect_x( b(3), -sku, sku);          %restrict skewness

[a0,a1,a2,a3,a4] = ak(m,s,T,S0,K,r);       %for all exercices prices simultaneously

%figure(1);
%plot(K,a0,K,a1,K,a2,K,a3,K,a4);

Cth= p0 * a0 + p3 * a3 + p4 * a4;

C=100*C;
Cth=100*Cth;

figure(1);
plot(K,Cth,K,C);