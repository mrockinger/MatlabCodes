% Example how to generate data when no explicit formula for CDF exists.
% Here illustrated in the case of a normal density

xa = -4;    % lower bound of support
dx = 0.01;  % increment
xb = +4;    % upper bound of support

x = xa:dx:xb;
fx = 1/sqrt(2*pi)*exp(-0.5*x.^2);  % here use normal density

plot(x,fx);
pause;

Fx = cumsum(fx)*dx; % the ad-hoc cdf

plot(x,normcdf(x,0,1),x,Fx); % verify that fit is good
pause;

NSim=1000;  % number of simulations
RandFx=zeros(NSim,1); % a placeholder for simulated data
u=rand(NSim,1); % simulate all numbers once and for all

for i=1:NSim
    y=Fx-u(i); 
    e=y>0;
    [m,idx]=max(e); %idx is first encounter where Fx>u(i)
    RandFx(i)=x(idx); % F^(-1)(x)
end;

hist(RandFx,10);



