% Random walk simulation
T=100;
Nsim=1;
x=rand(T,Nsim);
x=x-0.5
z=zeros(T,Nsim);
x=-(x<z)+(x>=z);
z1=zeros(1,Nsim)
x=[z1;x];
x=cumsum(x);
h=plot(x);
set(h,'LineWidth',2);
xlabel('Number of steps');
ylabel('S_N');;
title('101 steps of a random walk');
