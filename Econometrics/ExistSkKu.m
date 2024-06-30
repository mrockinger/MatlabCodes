function ExistSkKu(rt);
% Verifies the existence of skewness and kurtosis in vectorized format
clc;

load rtSP;
rt=rtSP;

[mi,id]=min(rt)

rt=[rt(1:4640); rt(4650:end)];

T=size(rt,1);
s=1:T; s=s';

mt1=cumsum(rt)./s;
mt2=cumsum(rt.^2)./s;
mt3=cumsum(rt.^3)./s;
mt4=cumsum(rt.^4)./s;

vt=mt2-mt1.^2; 
stdt=sqrt(vt);

sk=( mt3 - 3*mt2.*mt1 + 2*mt1.^3 )./( stdt.^3 );
ku=( mt4 - 4*mt3.*mt1 + 6*mt2.*mt1.^2 - 3*mt1.^4 )./( stdt.^4 );

subplot(2,1,1)
plot(sk);
title('skewness');

subplot(2,1,2)
plot(ku);
title('kurtosis');

