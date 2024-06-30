function test_ml();
% tests the maximum likelihood module
T=8100;
Rt=1+2*randn(T,1);

sum(lik([0.8;1.4],Rt))

beta0=[3; 1];
lb=[-10; 0.001];
ub=[+10; +100];


options=optimset('Diagnostics','on','Display','iter');

[beta,stderr2,vc,logl]=Max_lik(@lik,beta0,'Hessian',[],[],[],[],lb,ub,[],options,Rt);
[beta,stderr1,vc,logl]=Max_lik(@lik,beta0,'Sandwich',[],[],[],[],lb,ub,[],options,Rt);

[stderr1 stderr2]

function l=lik(beta,Rt);
% evaluates a column vector of log-likelihoods
mu =beta(1);
sig=beta(2);
l=-0.5*log(2*pi)-log(sig)-0.5*( (Rt-mu)/sig ).^2;
