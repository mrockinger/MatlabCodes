function plotCaFuNorm();
% plots the real part of the normal 
%characteristic function

t=-10:0.1:10;
mu=1;
sig=2;
b=CaFuNorm(t,mu,sig);
%plot(t,real(b));
%title('real part of ca norm');
%pause

size(b)

x=-10:0.1:10; x=x';
eitx=exp(-i*kron(t,x));

igd=matmul(eitx,b);
a=sum(igd,2)*(x(2)-x(1));
size(x)
size(real(a))
plot(x,real(a)')



function b=CaFuNorm(t,mu,sig)
% evaluates the characteristic function of a norma
b=exp(i*mu.*t-0.5*(sig.^2).*(t.^2));