function TestGMM(rt);
% GMM.m
% computes moments with the Generalized Method of Moments
% The GMM code comes from T. Roncalli.
clc;

load rtSP;
rt=rtSP;


beta0=[ mean(rt); std(rt) ; 0.1 ;0.1] ; %some initial values (mean, std, skewness, kurtosis)
disp('stating value ');
beta0

disp('verification that first two moments take value of 0');
m=Four_Mom(beta0,rt);
mean(m)
disp('Type any key to continue');
pause

lb=[-10; 0.0001; -100; 0.0001]; % lower bounds for moments
ub=[ 10; 10; 100; 100]; %upper bounds for moments

GMMLags=5;
GMMiter=1000;
GMMtol1=1e-5;
GMMtol2=1e-5;

options=optimset('Diagnostics','on','Display','iter');


[ beta, stderr, covbeta, Qmin, test, ptest ] = GMM(@Four_Mom,beta0,...
                              [],[],[],[],lb,ub,[],options,...
                              GMMLags, GMMiter, GMMtol1, GMMtol2,rt);
%    moments_generalises(start,rt,mmg_lags,mmg_iterations,mmg_tolerance1,mmg_tolerance2);


% the following is a Wald test for skewness and kurtosis being jointly nil

g=[ 0 0 1 0 ; 0 0 0 1];
sigma=g*covbeta*g';
ksi = (g*beta)'*inv(sigma)*(g*beta);
%ksi chi2cdf(ksi,2)


function mt = Four_Mom(beta,rt);
% here involves first four moments and first order autocorrelation. Returns moments 
% each observation. Hence result is a matrix 
% beta[1)=mean 
% beta(2)=standard deviation (not !! variance) 
% beta(3)=skewness 
% beta(4)=excess kurtosis

T=size(rt,1);
mt = zeros(T,4); 

mu  = beta(1);
sig = beta(2);
sk  = beta(3);
ku  = beta(4);

r2      = rt-mu;             % center data
mt(:,1) = r2;                % centered observations
mt(:,2) = (r2.^2)-sig^2;     % square sigma to get a variance
r3      = r2/sig;           % studentised data
mt(:,3) = (r3.^3) - sk;      % skewness
mt(:,4) = (r3.^4) - 3 - ku;  % excess kurtosis

 