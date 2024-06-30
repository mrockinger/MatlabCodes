function discretIllustrationBS()
% illustrates how discretization of call prices yields RND.
% also indicates limitations 
s=100;
r=0.05;
sigma=0.2;
k=70:2:130; k=k';
tau=60/365;
bsc = BSCall(s,k,sigma,tau,r);
figure(1)
plot(k,bsc)
title('Black-Scholes call option prices')
nk=rows(k);
dbsc=bsc(2:nk)-bsc(1:nk-1); dbsc/(k(2)-k(1));
nkm1=nk-1;
ddbsc=(dbsc(2:nkm1)-dbsc(1:nkm1-1))/(k(2)-k(1))^2;

eps=randn(nk,1);
bsc2=bsc + 0.01*mean(bsc)*eps;

nk=rows(k);
dbsc=bsc2(2:nk)-bsc2(1:nk-1); dbsc/(k(2)-k(1));
nkm1=nk-1;
ddbsc2=(dbsc(2:nkm1)-dbsc(1:nkm1-1))/(k(2)-k(1))^2;

figure(2)
h=plot(k(2:nkm1),ddbsc,k(2:nkm1),ddbsc2)
set(h,'LineWidth',2,{'LineStyle'},{'--';':'})
legend('RND with clean option prices','RND with perturbated option prices')
title('RND original and disturbed')

