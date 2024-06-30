function test_ml4(Rt);
% tests the maximum likelihood module
% estimates a GARCH(1,1) under the assumption that residuals are normal

%load rt3;
%Rt=rt3(2:end,:);

load rtSP;
Rt=rtSP;

m=mean(Rt);
Rt=Rt-m;

% parameters for GARCH: w, alpha, beta
beta0=[0.1; 0.05; 0.9];
lb=[0.001; 0.0001; 0.0001 ]; 
ub=[1; 0.3; 0.99];
A=[0 1 1]; %Ax<b
b=1;

options=optimset('Diagnostics','on','Display','iter');

[beta,stderr1,vc,logl]=Max_lik(@likGARCH11,beta0,'Hessian',A,b,[],[],lb,ub,[],options,Rt);
[beta,stderr2,vc,logl]=Max_lik(@likGARCH11,beta0,'Sandwich',A,b,[],[],lb,ub,[],options,Rt);

[stderr1 stderr2]

function l=likGARCH11(beta,rt);
% Computes the likelihood vector of a GARCH11
% presumes that Rt is a vector of centered returns 
w     =  beta(1); 
alpha =  beta(2); 
beta  =  beta(3);

rt2   = rt.^2;
[T,K] = size(rt2); 
ht    = zeros(T,1);
ht(1) = sum(rt2)/T;

for i=2:T
    ht(i)=w + alpha*rt2(i-1) + beta * ht(i-1);
end

sqrtht  = sqrt(ht);
x       = rt./sqrtht;

l = -0.5*log(2*pi) - log(sqrtht) - 0.5*(x.^2);

