% Simdtree1.prg
% simulates a tree with N steps for graphical purposes.
% Parameters are S0, r, sig, T, N

r=0.1;
r=log(1+r);
T=180/365;
N=20;
D=T/N;
sig=0.2;
S0=100;
NSim=N;

uN=exp(sig*sqrt(D));
dN=exp(-sig*sqrt(D));
dv=dN*ones(1,NSim);
uv=uN*ones(1,NSim);

A=zeros(NSim,NSim);
for i = 1:NSim;
    va=[dv(1:i) uv(i+1:NSim)]; 
    A(i,:)=va
end
A=A';
St1=S0*cumprod(A);

B=zeros(NSim,NSim);
for i = 1:NSim;
    vb=[uv(1:i) dv(i+1:NSim)]; 
    B(i,:)=vb
end
B=B';
St2=S0*cumprod(B);

St=100*ones(1,2*NSim);
St=[St;[St1 St2]];
plot(St);
%axis([0 2 80 130]);
set(findobj('Type','line'),'Color','k')