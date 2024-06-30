% plot2.m  again a plot
%
x=randn(10,1);
y=1+0.2*randn(10,1);
plot(x,y)
get(gcf,'Position') %gcf is generic graphics window handle
