function MixDRet()
% constructs the density for a mixed normal 
clc; cla; clear;

Param=[1 0.01 1 1; 1 1 1 1; 1 2 1 1; 1 2 1 10];
[nP,c]=size(Param);

ResM=[];

hold on
for i=1:4;
    mr=Param(i,1);
    sr=Param(i,2);
    mn=Param(i,3);
    sn=Param(i,4);

if i==1    
   ng=(0.0001:0.01:25);
elseif i==2
   ng=(0.0001:0.1:25);   
elseif i==3
   ng=(0.0001:0.1:25);   
else
   ng=[(0.0001:0.0001:2) (2:0.2:50)];   
end

fn = 1./ ( sqrt(2*pi).*sn.*ng ) .* exp( -0.5*( ( log(ng)-mn )./sn ).^2 );
%plot(ng,fn)
%axis([0 0.2 0 50])
%return

% now go matrix
if i==1
    r=(-1:0.01:25)';
elseif i==2
    r=(-4:0.1:35)';
elseif i==3
    r=(-12:0.1:35)';
else
    r=(-3.5:0.01:3.5)';
end
rM = kron(r,ones(size(ng)));
nM = kron(ng,ones(size(r)));
grn= 1./ ( sqrt(2*pi).*sr.*sqrt(nM) ) .* exp( -0.5*( ( rM - mr.*nM ).^2 )./( (sr.^2).*nM ) );
fnM= kron(fn,ones(size(r)));
fr=sum( grn.*fnM,2)*( ng(2)-ng(1) );

plot(r,fr)

s=sum(fr)*(r(2)-r(1))
fr=fr/s; % scale so that probability mass equals 1

% compute selected moments.
m1 =sum(r.*fr)*(r(2)-r(1))
v  =sum( ((r-m1).^2) .*fr)*(r(2)-r(1))
se=sqrt(v);
sk =sum( ((r-m1)./se).^3 .*fr)*(r(2)-r(1))
ku =sum( ((r-m1)./se).^4 .*fr)*(r(2)-r(1))
Resv=[s m1 v se sk ku]';
ResM=[ResM Resv];

end
hold off

for i=1:6
    fprintf('%8.2f %8.2f %8.2f %8.2f', ResM(i,:))
end
