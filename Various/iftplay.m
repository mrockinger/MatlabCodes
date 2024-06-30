function iftplay()
% iftplay.m implents the inverste Fast Fourier Transform
clc;
u=-6:0.1:6; u=u';

mu=4;
sig=1;

figure(1)
f=cf(u,mu,sig);

plot(u,real(f))
title('Real part of characteristic function. N(4,1)')

xL=0;
xU=10;
N=2^10;
Dx=(xU-xL)/N;
xk=xL:Dx:xU;

Du=(2*pi)/(N*Dx);
uL= 0;
uU= Du*N;
u = uL:Du:uU;

cf1  = cf(u,mu,sig);
phi = fft( exp(-i*Du*(N-1)*xL) .* cf1 );
phi = real(phi)/pi;



%plot(xk,phi)
%title('Discrete Inverse FFT of the normal cf')


% exact integration at one point and downward adjustment
uL = 0;
uU = 10;
Ne  = 1000;
du = (uU-uL)/Ne;
u  = uL:du:uU;

cf2    = cf(u,mu,sig);
truef0 = sum(real( exp(-i*u*xL) .*cf2 ))*du/pi;
df     = phi(1)-truef0;

phifin = phi-df;

figure(2)
h=plot(xk,phi,xk,phifin)
set(h,'LineWidth',2,{'LineStyle'},{'--';'-.'})
set(h,{'Color'},{'m';'m'})
xlabel('x')
ylabel('f(x)')
legend(h,'Density out of FFT','Corrected density')
title('Discrete Inverse FFT of the normal cf')
axis([0 10 0 1])

dnorm=1/(2*pi)*exp(-0.5*((xk-mu).^2)/(sig^2));
diff=dnorm-phifin;
max(diff)






function f=cf(u,mu,sig)
% implementation of characteristic function
% here implements the Normal mu , sig^2
f=exp( i.*mu.*u - 0.5*(sig^2).*(u.^2) );
