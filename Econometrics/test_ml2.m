function test_ml2();
% tests the maximum likelihood module
% here uses interest rate data and estimates the mean and variance

%load rt3;
%Rt=rt3(2:end,:);

load rtSP;
Rt=rtSP;

beta0=[-5; 1];
lb=[-10; 0.001];
ub=[+10; +100];

options=optimset('Diagnostics','on','Display','iter');

[beta,stderr1,vc,logl]=Max_lik(@lik,beta0,'Hessian',[],[],[],[],lb,ub,[],options,Rt);
[beta,stderr2,vc,logl]=Max_lik(@lik,beta0,'Sandwich',[],[],[],[],lb,ub,[],options,Rt);


[stderr1 stderr2]

function l=lik(beta,Rt);
% evaluates a column vector of log-likelihoods
mu =beta(1);
sig=beta(2);
l=-0.5*log(2*pi)-log(sig)-0.5*( (Rt-mu)/sig ).^2;
