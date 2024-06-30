function TraceLognormal()
% trace lognormal
% traces a lognormal 

Resfn=[]; %store info for later on
Resng=[];

Param=[1 0.01 1 1; 1 1 1 1; 1 2 1 1; 1 2 1 10];
[nP,c]=size(Param);

hold on
for i=1:nP;
    mr=Param(i,1);
    sr=Param(i,2);
    mn=Param(i,3);
    sn=Param(i,4);

if i<=3    
   ng=(0.0001:0.1:10);
else
   ng=(0.0001:0.1:10);   
end

fn = 1./ ( sqrt(2*pi).*sn.*ng ) .* exp( -0.5*( ( log(ng)-mn )./sn ).^2 );

plot(ng,fn)
axis([0 10 0 0.5])
end
hold off
return

Resfn=[Resfn; fn];
Resng=[Resng; ng];

% now go matrix
r=(-6:0.1:25)';
rM = kron(r,ones(size(ng)));
nM = kron(ng,ones(size(r)));
grn= 1./ ( sqrt(2*pi).*sr.*sqrt(nM) ) .* exp( -0.5*( ( rM - mr.*nM ).^2 )./( (sr^2).*nM ) );
fnM= kron(fn,ones(size(r)));
fr=sum( grn.*fnM,2)*(ng(2)-ng(1));
figure(i)
plot(r,fr)

sum(fr)*(r(2)-r(1))
end

figure
hold on
for i=1:4
plot(Resng(i,:),Resfn(i,:))
end