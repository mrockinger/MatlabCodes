function xst=Studentize(x)
% ensures that each column of x has mean 0 and variance 1
m=mean(x)';
s=std(x)';
T=rows(x);
ov=ones(T,1);
xst= (x-kron(m,ov))./kron(s,ov);
