function BasicGraph()
% do not confuse with basic instincts
echo on
u=randn(5,2) % 5 uniform random numbers in two columns
u(:,2)=u(:,1)+0.3*u(:,2)

plot(u(:,1),u(:,2),'s',...
 'MarkerEdgeColor','k',...
 'MarkerFaceColor','g',...
 'MarkerSize',10)
title('scatterplot of u(:,1) against u(:,2)');
% in the plot instruction 's' means symbols="Marker"
[ma,ma_idx]=max(u)

max(u) % notice that Matlab functions are overloaded

[mi,mi_idx]=min(u)

mu=mean(u)

sum(u)/size(u,1)

median(u)

sort(u)

sortrows(u)

uu = u - kron( mu, ones(size(u,1),1)); %this is truly stupid
mean( uu.^2 )

T=size(u,1);
mean(uu.^2)*(T-1)/T

v=var(u)

sqrt(v)

std(u)

mu=0.1;
sig=0.2;
S0=100;
x=mu+sig*u;
x1=exp(x);
x1=[ones(1,2); x1]; % first observation
St1=100*cumprod(x1)

x=[zeros(1,2); x]; % again for first observation
St2=100*exp(cumsum(x))

A=[1 2 3; 2 2 4; 2 1 2]

prod(A)

cumprod(A)

dx=0.1

x=-5:dx:5; % some fine grid
y=1/sqrt(2*pi)*exp(-0.5*x.^2);
figure(1)
subplot(2,1,1)
plot(x,y);
cdf=cumtrapz(y)*dx;
subplot(2,1,2)
plot(x,cdf);
trapz(y)*dx
  