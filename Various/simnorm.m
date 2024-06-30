% simulation of a normal sample
N=1000;
x=randn(N,1);

clf;

subplot(2,1,1);
plot(x,'o');
xlabel('x');
ylabel('normal draws');
%axis([0 1000 -4 4]);
set(findobj('Type','line'),'Color','k')

subplot(2,1,2);
hist(x,20);
xlabel('x');
ylabel('histogram');
