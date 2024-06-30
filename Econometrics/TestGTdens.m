function TestGTdens()
% constructs the Generalized Student t over a given support
% then computes quantiles and numerical expected shortfall

dz=0.001;
z=-6:dz:6; z=z';
lda=0;
eta=30;
g1 = GTdens(z,lda,eta);

lda=0;
eta=5;
g2 = GTdens(z,lda,eta);

lda=-0.8;
eta=5;
g3 = GTdens(z,lda,eta);



figure(1)
plot(z,g1,z,g2,z,g3)
title('Density of Skewed Student t')
text(0,0,'\lambda=0, \eta=30')
text(0,0.1,'\lambda=0, \eta=5')
text(0,0.2,'\lambda=-0.8, \eta=5')
cg1=cumsum(g1);
figure(2)
plot(z,cg1*dz)
title('CDF of asymmetric Student t')

% computes numerically the VaR and Expected shortfall
size(cg)
conflevel=0.01;
conflevel.*ones(size(cg))-cg*dz<0
[level idx]=min((conflevel.*ones(size(cg))-cg*dz)>0);
idx
cg(idx-10:idx+10)*dz
%compute expected shortfall
sum( z(1:idx).*cg(1:idx)*dz )/sum( cg(1:idx)*dz )

