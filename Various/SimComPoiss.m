function SimCompPoiss(lambda,T,jumpType,m,v)
% returns a trajectory of a compound Poisson process
% the returned trajectories are with drift and compensated, i.e. martingale behavior
% input:
% lambda the jump intensity
% T the length of time
% jumpType is the type of required jumps.
%    if equals to 0 then use normal distribution with mean m and variance v
%    if equals to 1 
% 
% general parameters
% N the number of discrete steps over which to decompose [0,T]



lambda=10; % number of jumps per time interval 
T=1; % time interval, here taken to be 1 (year or whatever)
m=5; % mean of normal jump size
s=10; % standard deviation of jumps

mu=-lambda*m; % compensator

nJ  = SimPoiss(lambda*T,1)       % one draw to determine number of jumps
tv = 0:0.001:T; tv=tv';          % time
N=rows(tv);
jv=1:1:N;                        % discrete grid
jT = unifUnique(nJ,1,N);         % nJ unique uniformly distributed numbers between 1 and N


u   = rand(nJ,1); dur=-log(u)/(lambda*T); % draw the durations between jumps (exponential durations)
Tj  = cumsum(dur) % times when a jumps is to occur
Js  = m+randn(nJ,1)*s
Js  = cumsum(Js);


Pp = mu*tv;

% convert time into numbers 1...N that are useful for the display
N=rows(tv);
j=1:1:N;
Tj=floor( (N-1)/T*Tj + 1 ); %times of jumps as an index
Xt=zeros(Tj,1);
for i=2:nJ
    Tj(i)-Tj(i-1)+1
    Xt=[Xt;ones(Tj(i)-Tj(i-1)+1,1)*Js(i-1)];    
end

Pp(Tj(i))=Pp(Tj(i))-Js(i)
plot(tv,Xt)




