function [tv,Xt,CXt]=SimCompPoiss(lambda,T,jumpType,m,v)
% returns a trajectory of a compound Poisson process
% the returned trajectories are with drift and compensated, i.e. martingale behavior
% input:
% lambda the jump intensity
% T the length of time
% jumpType is the type of required jumps.
%    if equals to 0 then use normal distribution with mean m and variance v
%    if different from 0 then use constant jumps, of size jumpType 
% 
% general parameters
% N the number of discrete steps over which to decompose [0,T]
N=1000;

if jumpType == 0; %normal case
    mu=lambda*m; % compensator
else
    mu=lambda*jumpType;
end

nJ  = poissrnd(lambda*T,1,1)       % one draw to determine number of jumps
dt=T/N;
tv = 0:dt:T; tv=tv';          % time
N=rows(tv);
jv=1:1:N; jv=jv';                        % discrete grid
u   = rand(nJ,1); 
Xt  = zeros(N,1);

if nJ>0
    jIdx = unifUnique(nJ,1,N);         % nJ unique uniformly distributed numbers between 1 and N
    if jumpType==0
        jS = m + randn(nJ,1)*sqrt(v); % the jump size
    else
        jS = jumpType*ones(nJ,1);
    end
    
    for i=1:nJ
        Xt(jIdx(i):N) = Xt(jIdx(i):N) + jS(i)*ones(N-jIdx(i)+1,1);
    end
    
end

CXt = Xt - mu*tv;
