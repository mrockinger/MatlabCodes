% Stat1.m
%
% Simulates draws from a normal density and an exponential and 
% verify that the mean satisfies the consistency property

NSim=100;
lda=0.5;
mu=2;

% draw from exponential
u=rand(NSim,1);
expo=-1/lda*log(u);

% draw from normal
norm=mu+randn(NSim,1);

s=1:NSim; s=s';
norm=cumsum(norm)./s;
expo=cumsum(expo)./s;

plot(s,norm,'-d',s,expo,'-p');

% simulates a sample of size Msim of an exponential distribution
