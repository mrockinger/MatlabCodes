% StatCompDeltaSim.prg
% simulates geometric brownian motions and traces them as 
% an exercise of comparative static
% Here we simulate a set of brownian motions simultaneously.
% Parameters are S0, r, sig, T, N, K.

r=0.1;
r=log(1+r);
T=100/365;
N=1000;
D=T/N;
sig=0.2;
S01=100;
S02=150;
NSim=10;
K=100;


eps=randn(N,NSim);

% comparative statics concerning S_0
Xt = (r-0.5.*(sig^2)).*D + sig .* sqrt(D).* eps;
Xt=[zeros(1,NSim); Xt(2:N,:)];
Xt=cumsum(Xt);

idx=1:N;
St1=S01*exp(Xt);
St2=S02*exp(Xt);

% comparative statics concerning volatility
sig1=0.2;
sig2=0.5;
Xt1 = (r-0.5.*(sig1^2)).*D + sig1 .* sqrt(D).* eps;
Xt1=[zeros(1,NSim); Xt1(2:N,:)];
Xt1=cumsum(Xt1);

Xt2 = (r-0.5.*(sig2^2)).*D + sig2 .* sqrt(D).* eps;
Xt2=[zeros(1,NSim); Xt2(2:N,:)];
Xt2=cumsum(Xt2);
idx=1:N;
St3=S01*exp(Xt1);
St4=S01*exp(Xt2);


subplot(2,2,1);
plot(idx,St1,[N;N],[K;0]);
set(findobj('Type','line'),'Color','k')
axis([0 1010 0 200]);
title('Case of S_0=100');

subplot(2,2,2);
plot(idx,St2,[N; N],[K;0]);
set(findobj('Type','line'),'Color','k')
axis([0 1010 0 200]);
title('Case of S_0=150');

subplot(2,2,3);
plot(idx,St3,[N;N],[K;0]);
set(findobj('Type','line'),'Color','k')
axis([0 1010 0 200]);
title('Case of \sigma=0.2');

subplot(2,2,4);
plot(idx,St4,[N; N],[K;0]);
set(findobj('Type','line'),'Color','k')
axis([0 1010 0 200]);
title('Case of \sigma=0.5');