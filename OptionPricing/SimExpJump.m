function [tv,Jv]=SimExpJump(N,T,lda,eta)
% Description:
%   Creates a grid for time (potentially useful for a plot). Also creates
%   the associated jump vector. This jump vector is actually a step function.
% Input:
%   N the number of time increments
%   T the total time interval [0, T]
%   lda is the jump intensity
%   eta is the jump size (assumed to be an exponential)
% Output:
%  tv a 
%  Jv the step function containing the jumps

% create the time vector
dt = T/(N-1);
tv = (0:dt:T)';                % time. This is useful if a plot is required
jv = (1:N)';                   % discrete grid to contain jump indices
nJ = poissrnd(lda*T,1,1);   % one draw to determine number of jumps
%N
%size(tv)

u   = rand(nJ,1); 
Jv  = zeros(N,1);

if nJ>0
    jIdx = unifUnique(nJ,1,N);         % nJ unique uniformly distributed numbers between 1 and N
    jS   = exprnd(eta,nJ,1);           % the jump size is exponential
    
    for i=1:nJ
        Jv(jIdx(i):end) = Jv(jIdx(i):end) + jS(i)*ones(N-jIdx(i)+1,1);
    end
    
end
