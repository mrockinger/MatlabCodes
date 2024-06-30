function LevyProcess1()
% now combine a compensated Poisson process with an arithmetic BM

%Parameters for the compound Poisson process
lambda=3;
T=10; 
jumpType=0; % corresponds to N(m,v) jumps
%jumpType=1; % jumps with constant size (Size=jumpType)
m=3;
v=9;

% Parameters for the arithmetic BM
mu=2;
sig=4;

[Tv,Xt,CXt]=SimCompPoiss(lambda,T,jumpType,m,v);

nT  = rows(Tv);
dT  = Tv(2)-Tv(1);
dZt = mu*dT + sig*sqrt(dT)*randn(nT-1,1); %increments of BM
Zt  = cumsum(dZt); %sum of increments yields BM
Zt  = [0;Zt];

if jumpType==0;

figure(1)
plot(Tv,Xt)
title('Compound Poisson process')

size(Tv)
size(CXt)
size(Zt)
figure(2)
plot(Tv,CXt+Zt)
title('Aritmetic BM + compensated compound Poisson process')

else

figure(1)
plot(Tv,Xt)
title('Poisson process with jump size of 1')

figure(2)
plot(Tv,CXt)
title('Compensated Poisson process with jump size of 1')    
end
