function RandFx=SimPoiss(lambda,NSim)
% simulation of a set of N Poission distributed random variates with
% parameter lambda
M=20; % simulate Poisson rv based on a discretization of the CDF involving M steps
kv=0:1:M; kv=kv';
sP = -lambda + kv*log(lambda) - gammaln(kv+1); 
sP = exp(sP);
FxP=cumsum(sP);
FxP=[FxP;1];
plot(FxP)

% now imbed problem into general inverse CDF simulator
RandFx=zeros(NSim,1); % a placeholder for simulated data
u=rand(NSim,1); % simulate all numbers once and for all

for i=1:NSim
    y=FxP-u(i); 
    e=y>0;
    [m,idx]=max(e); %idx is first encounter where Fx>u(i)
    
    RandFx(i)=kv(idx); % F^(-1)(x) gives the k we are seeking
end;