% cdf and inverse cdf 
x=-4:0.1:4;
dens=1./sqrt(2.*pi).*exp(-0.5*x.^2);
cdf=cumsum(dens)*(x(2)-x(1));

N1l=25;
x1l=x(1:N1l);
y1l=ones(N1l,1)*cdf(N1l);
N2l=N1l+2;
x2l=x(1:N2l);
y2l=ones(N2l,1)*cdf(N2l);
x3l=x(N1l)*ones(2,1);
y3l=[cdf(N1l);0];
x4l=x(N2l)*ones(2,1);
y4l=[cdf(N2l);0];

N1u=40;
x1u=x(1:N1u);
y1u=ones(N1u,1)*cdf(N1u);
N2u=N1u+2;
x2u=x(1:N2u);
y2u=ones(N2u,1)*cdf(N2u);
x3u=x(N1u)*ones(2,1);
y3u=[cdf(N1u);0];
x4u=x(N2u)*ones(2,1);
y4u=[cdf(N2u);0];


subplot(2,1,1); 
title('Inverse cdf method');
plot(x,cdf,x1l,y1l,x2l,y2l,x3l,y3l,x4l,y4l,x1u,y1u,x2u,y2u,x3u,y3u,x4u,y4u);
set(findobj('Type','line'),'Color','k')
xlabel('x');
ylabel('cdf');

subplot(2,1,2); 
y3l(1)=0.4;
y4l(1)=0.4;
y3u(1)=0.4;
y4u(1)=0.4;
plot(x,dens,x3l,y3l,x4l,y4l,x3u,y3u,x4u,y4u);
xlabel('x');
ylabel('density');
set(findobj('Type','line'),'Color','k')