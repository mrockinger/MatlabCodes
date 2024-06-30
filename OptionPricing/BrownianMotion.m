% simulates brownian motion
T=10000;
x=randn(T,2);
x=[0 0;x];
x=cumsum(x);
h=plot(x(:,1),x(:,2));
set(h,'LineWidth',1.5);
