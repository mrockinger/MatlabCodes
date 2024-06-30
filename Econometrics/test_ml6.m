function test_ml6(Rt);
% tests the maximum likelihood module
% estimates a GARCH(1,1) assuming that residuals are Student-t
% incorporates asymetry for good and bad news

load rtSP;
Rt=rtSP;

m=mean(Rt);
Rt=Rt-m;

% parameters for GARCH: w, alpha, beta, nu
beta0=[0.01; 0.1; 0.1; 0.8; 4];
lb=[0; 0.0001; 0.0001; 0.0001; 2 ]; 
ub=[1; 0.3;    0.3;    0.99;   30];
A=[0 0.5 0.5 1 0]; %Ax<b
b=1;

[ lb beta0 ub]

options=optimset('Diagnostics','on','Display','iter');

[beta,stderr1,vc,logl]=Max_lik(@likGARCH11T,beta0,'Hessian',A,b,[],[],lb,ub,[],options,Rt);
[beta,stderr2,vc,logl]=Max_lik(@likGARCH11T,beta0,'Sandwich',A,b,[],[],lb,ub,[],options,Rt);

[stderr1 stderr2]

function l=likGARCH11T(b,rt);
% Computes the likelihood vector of a GARCH11
% presumes that Rt is a vector of centered returns 

w      = b(1); 
alphap = b(2); 
alphan = b(3);
beta   = b(4);
nu     = b(5);

rt2   = rt.^2;
ep    = rt>0;
en    = rt<=0;

[T,K] = size(rt2); 
ht    = zeros(T,1);
ht(1) = sum(rt2)/T;

x = w + alphap*rt2.*ep + alphan*rt2.*en;

for i=2:T
    ht(i) =  x(i-1) + beta * ht(i-1);
end

sqrtht  = sqrt(ht);
x       = rt./sqrtht;

%*** Plug in Student-t

np1= (nu+1)/2;
nh = nu/2;
l1 = gammaln(np1);
l2 = gammaln(nh);
l3 = log(sqrtht);
l4 = 0.5*log(nu*pi);
l5 = np1*log(1 + (1/nu) * (x.^2));

l=l1-l2-l3-l4-l5;
