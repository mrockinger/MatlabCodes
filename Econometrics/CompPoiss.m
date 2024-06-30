function CompPoiss()
% simulates a compound poisson process
% simplest version
% formula for process comes from...
lda= 1;
T  = 100;
jt = 0;
TT = 0;
Nj = 0;
while TT<T
        el=2*rndexp(lda,1)
        TT=TT+el;
        jt=[jt;TT];
        Nj=Nj+1;
end
%jt=jt(1:end-1);

% actual jumps now
Xt=[0; randn(Nj,1)];
Xt=cumsum(Xt);
stairs(jt,Xt)



function x=rndexp(lda,n)
x=-log(rand(n,1))/lda;