function TestSimCompPoiss()

lambda=30;
T=10; 
jumpType=0; % corresponds to N(m,v) jumps
%jumpType=1; % jumps with constant size (Size=jumpType)
m=3;
v=9;

[Tv,Xt,CXt]=SimCompPoiss(lambda,T,jumpType,m,v);

if jumpType==0;

figure(1)
plot(Tv,Xt)
title('Compound Poisson process')

figure(2)
plot(Tv,CXt)
title('Compensated compound Poisson process')

else

figure(1)
plot(Tv,Xt)
title('Poisson process with jump size of 1')

figure(2)
plot(Tv,CXt)
title('Compensated Poisson process with jump size of 1')    
end
