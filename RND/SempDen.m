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
% To run program, must have sk.fmt containing the skewness-kurtosi:s boundary
%
clc; 

%**********************************************************************
%*                   load data and set global parameters              *
%**********************************************************************
%globals for skewness+kurtosis boundaries
eps=0.0001; % distance from boundary
s6=sqrt(6);
s24=sqrt(24);

[fp, msg] = fopen('SK.txt','r');
if fp == -1
    disp(msg)
end
[sk,c]=fscanf(fp,['%f %f'],[2,8001]);
sk=sk';


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

eps=1e-07;
ml_sel=1;

titre='Estimation of semi-parametric risk neutral density for 250497';

param=[0.002   0.030   0.559   0.695 ;
       0.001   0.031   0.508   0.843 ;
       0.002   0.031   0.487   0.938 ;
       0.001   0.032   0.432   1.047 ;
       0.002   0.032   0.401   1.050 ;
       0.002   0.032   0.436   1.055 ];

ResPara=[];
sempm=[];

for selmat=1:6;

param0 = param(selmat,:)';
param0=param0+0.1; 

[ S0, C, K, r, rs, T] = initdata(selmat);
C=flipud(C);
K=flipud(K);

r  = log(1+r);   %continuously compounded interest rates for France
rs = log(1+rs);  %continuously compounded interest rates for Germany

%To make further computations easier
C = C*exp(r*T);

A = [];
b = [];
Aeq = [];
beq = [];
nonlcon = [];
options = optimset('Diagnostics','on','Display','iter');
lb = [-0.1; 0.025; 0.0001; 0.0001 ]; 
ub = [ 0.1; 0.04;    30;   30];

How_is_fit(param0,kumin,kumax,sk,ai,bi,T,S0,K,C)
%pause

[b,fmin,retcode,output,lambda,gr,hessian] =...
      fmincon(@obj,param0,A,b,Aeq,beq,lb,ub,nonlcon,options,...
                kumin,kumax,sk,ai,bi,T,S0,K,C);

fprintf('Result of optmum');
fprintf('maturity %8.4f', selmat);
b

% verify fit
How_is_fit(b,kumin,kumax,sk,ai,bi,T,S0,K,C)
ResPara=[ResPara; b'];
ResPara

b1=b;
b1(4)=rect_x(b1(4),kumin,kumax);       %map kurtosis into allowed domain%
sku=interk(b1(4),sk,ai,bi);                   %get upper bound for skewness%
b1(3)=rect_x(b1(3),-sku,sku);          %restrict skewness%
b1';

m  = b(1) ;
s  = b(2);
b0 = 1;
b1 = 0;
b2 = 0;
b3 = b(3);
b4 = b(4);


%recover the traditional risk neutral density for S_T*/
z    = 3.2:0.005:3.6; z=z';
mu   = log(S0) + ( m - 0.5 * s^2 ) * T;
si   = s*sqrt(T);
tmp  = ( log(z) - mu )/si;
semp = 1/si * ri_ne_de(tmp,[b0;b1;b2;b3;b4],kumin,kumax,sk,ai,bi);
semp = semp./z;

sempm = [sempm semp];

end


%_plegstr="30 days\00060 days\00090 days\000180 days\000270 days\0001 year";
%_plegctl={2 5 6 4};
%_pdate="";
%_plwidth=2;
%xtics(3.2,3.60,0.04,1);
%ytics(0,25,2,1);
figure(2)
size(sempm)
plot(z,sempm(:,1),z,sempm(:,2),z,sempm(:,3),z,sempm(:,4),z,sempm(:,5),z,sempm(:,6));

title('Estimated semi-parametric density function');
%xlabel('S]T[');
%ylabel('Density');

disp('the parameters are ');
ResPara

%******************************************************************************

function [ai,bi]=init_ab(sk,eps); %initializes parameters for linear interpolation
N  = rows(sk);
k  = sk(:,1);
s  = sk(:,2) - eps;

dk = k(2:N) - k(1:N-1);
ai = ( s(1:N-1) .* k(2:N) - k(1:N-1) .* s(2:N) ) ./ dk;
bi = ( s(2:N) - s(1:N-1) ) ./ dk;

%-----------------------------------------------------------------------------

function f = obj(b,kumin,kumax,sk,ai,bi,T,S0,K,C);

m  = b(1); 
s  = b(2);
p0 = 1; 
p3 = b(3); 
p4 = b(4);

p4  = rect_x( b(4), kumin, kumax);       %map kurtosis into allowed domain
sku = interk(p4,sk,ai,bi);               %get upper bound for skewness
p3  = rect_x( b(3), -sku, sku);          %restrict skewness

[a0,a1,a2,a3,a4] = ak(m,s,T,S0,K);       %for all exercices prices simultaneously

%plot(K,a0,K,a1,K,a2,K,a3,K,a4);
%pause

Cth= p0 * a0 + p3 * a3 + p4 * a4;

%plot(K,Cth,K,C);
%pause

dist = C - Cth;
f    = dist' * dist;
f = 10000*f

%-----------------------------------------------------------------------------

function [a0,a1,a2,a3,a4]=ak(m,s,T,S0,K);
%returns a0, a1, a2, a3, a4

st = s*sqrt(T);
sp = 1/sqrt(2*pi);
d1 = ( log(S0./K) + ( m + 0.5*(s.^2)).*T )./st;
d2 = d1 - st; 

f   = S0*exp(m.*T); 
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

%proc (1)=invrect(y,a,b);
%%x is a vector with elements ranging between ]a,b[;
%This procedure computes the inverse to get an x between -infty,+infty %
%local u,x;
%u= (y-a)/(b-a);
%x=ln( u/(1-u) );
%retp(x);
%endp;


%-----------------------------------------------------------------------------


function [S0,resC,resK,rFr,rDe,dte]=initdata(mat);
S0=3.3735;
resC=[0.0020    0.0030    0.0037    0.0053    0.0065    0.0076;
   0.0030    0.0044    0.0054    0.0075    0.0092    0.0108;
   0.0039    0.0056    0.0068    0.0093    0.0114    0.0135;
   0.0057    0.0081    0.0099    0.0135    0.0165    0.0195;
   0.0076    0.0105    0.0126    0.0173    0.0213    0.0248;
   0.0100    0.0139    0.0167    0.0224    0.0275    0.0321;
   0.0128    0.0180    0.0214    0.0284    0.0357    0.0417;
   0.0160    0.0222    0.0266    0.0349    0.0440    0.0515;
   0.0209    0.0290    0.0348    0.0456    0.0578    0.0680;
   0.0242    0.0336    0.0403    0.0523    0.0666    0.0786;
   0.0296    0.0402    0.0483    0.0631    0.0809    0.0964 ];

resK = [   3.4292    3.4568    3.4800    3.5286    3.5652    3.5998;
   3.4145    3.4344    3.4491    3.4813    3.5066    3.5306;
   3.4037    3.4175    3.4283    3.4502    3.4683    3.4853;
   3.3899    3.3973    3.4035    3.4158    3.4259    3.4351;
   3.3807    3.3841    3.3870    3.3936    3.3985    3.4022;
   3.3740    3.3747    3.3757    3.3779    3.3791    3.3796;
   3.3679    3.3662    3.3655    3.3643    3.3620    3.3594;
   3.3621    3.3583    3.3560    3.3519    3.3462    3.3409;
   3.3554    3.3489    3.3447    3.3370    3.3271    3.3181;
   3.3513    3.3433    3.3379    3.3285    3.3160    3.3046;
   3.3453    3.3357    3.3287    3.3162    3.2996    3.284];

rFr = [0.0333 0.0339 0.0345 0.0353 0.0357 0.0360]; rFr=rFr';
rDe = [0.0316 0.0316 0.0316 0.0322 0.0328 0.0334]; rDe=rDe';

dte = [0.0822 0.1644 0.2466 0.4932 0.7397 1.0000]; dte=dte';

resC = resC(:,mat),
resK = resK(:,mat),
rFr  = rFr(mat),
rDe  = rDe(mat),
dte  = dte(mat);

%-----------------------------------------------------------------------------

function How_is_fit(b,kumin,kumax,sk,ai,bi,T,S0,K,C)
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

[a0,a1,a2,a3,a4] = ak(m,s,T,S0,K);       %for all exercices prices simultaneously

%figure(1);
%plot(K,a0,K,a1,K,a2,K,a3,K,a4);

Cth= p0 * a0 + p3 * a3 + p4 * a4;

C=100*C;
Cth=100*Cth;

disp('Actual and new fit of options');
[K C Cth]
disp('now plot this shit');

Pb=[K C K Cth];
n=rows(Pb);
for i=1:n
fprintf('Pb %15.6f %15.6f %15.6f %15.6f\n',Pb(i,:));
end

%figure(1);
plot(K,Cth,K,C);