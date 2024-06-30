function test_ml3(Rt);
% tests the maximum likelihood module
% here uses SP500 returns and estimates the mean and variance
% estimates the parameters for a Student-t

 %Rt=2*trnd(5,8105,1);
 %disp('theoretical mean and variance and empirical counterpart');
 %[[mean(Rt); var(Rt)] [0;4*(5/(5-2))] ]


load rtSP;
Rt=rtSP;

beta0=[0; 1; 7]; %new parameter is nu
lb=[-10; 0.001; 2]; % and nu should be larger than 2
ub=[+10; +100; 100];

options=optimset('Diagnostics','on','Display','iter');

[beta,stderr1,vc,logl]=Max_lik(@likStud,beta0,'Hessian',[],[],[],[],lb,ub,[],options,Rt);
[beta,stderr2,vc,logl]=Max_lik(@likStud,beta0,'Sandwich',[],[],[],[],lb,ub,[],options,Rt);

[stderr1 stderr2]

function l=likStud(beta,Rt);
% evaluates a column vector of log-likelihoods
mu = beta(1);
sig= beta(2);
nu = beta(3);

np1= (nu+1)/2;
nh = nu/2;
l1 = gammaln(np1);
l2 = gammaln(nh);
l3 = log(sig);
l4 = 0.5*log(nu*pi);
Rtc= (Rt-mu)/sig;
l5 = np1*log(1 + (1/nu) * (Rtc.^2));

l=l1-l2-l3-l4-l5;
