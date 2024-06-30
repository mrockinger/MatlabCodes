% TestCallSim
% tests if the module that computes the call price price 
% via simulation works.
S0=100;
K=100;
sig=0.2;
r=0.1;
r=log(1+r);
T=100/365;
N=1000;
MSim=100000; %Number of paths to construct option price

tic;
CallSim(S0,K,sig,r,T,N,MSim);
t=toc;
t