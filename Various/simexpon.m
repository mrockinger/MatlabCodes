% simulation of an exponentially distributed sample
N=500;
lambda=0.5
u=rand(N,1);
x=-log(u)/lambda;

subplot(2,1,1);
plot(x,'o');
xlabel('x');
ylabel('exponential draws');
axis([0 500 0 15]);
set(findobj('Type','line'),'Color','k')

subplot(2,1,2);
hist(x,10);
xlabel('x');
ylabel('histogram');
