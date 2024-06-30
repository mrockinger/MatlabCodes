function test_mlmix1();
% tests the maximum likelihood module
% estimates the parameters for a mixture of distributions
% not the most stable code


load rtSP;
Rt=rtSP;

beta0=[0; 1; -1; 2; 0.1]; % mu1 sig1 mu2 sig2 a
beta0(1)=mean(Rt);
beta0(3)=mean(Rt);
beta0(2)=std(Rt);
beta0(4)=2*std(Rt);

lb=[-10; 0.001; -10; 0.001;0.0001]; 
ub=[+10; +100; 10; +100; 0.9999];

options=optimset('Diagnostics','on','Display','iter');
A=[0 1 0 -1 0];
b=0;

[beta,stderr1,vc,logl]=Max_lik(@likMix,beta0,'Hessian',A,b,[],[],lb,ub,[],options,Rt);
[beta,stderr2,vc,logl]=Max_lik(@likMix,beta0,'Sandwich',A,b,[],[],lb,ub,[],options,Rt);

[stderr1 stderr2]

function l=likMix(beta,Rt);
% evaluates a column vector of log-likelihoods
mu1 = beta(1);
sig1= beta(2);
mu2 = beta(3);
sig2= beta(4);
a   = beta(5);

Rtc1= (Rt-mu1)/sig1;
Rtc2= (Rt-mu2)/sig2;

f1=1/(sqrt(2*pi)*sig1)*exp(-0.5*(Rtc1.^2));
f2=1/(sqrt(2*pi)*sig2)*exp(-0.5*(Rtc2.^2));

l=log( a.*f1 + (1-a).*f2 );
