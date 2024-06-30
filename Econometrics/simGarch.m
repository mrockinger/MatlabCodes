function simGarch()
% simulates for a given horizon trajectories and fits a smoothed density
% the density is a kernel fit
% parameters correspond to the last 3000 observations of the FTSE 1000
clc;
cla;
m =0.0188;
w =0.0089;
ap=0.0047;
an=0.0854;
be=0.9369;
nu=14.8537;
sig0 = 0.9557;
S0   = 4357.53;

NbD  = 20; %time horizon in days
NSim = 10000;

xt = trnd(nu,NbD,NSim); %simulate following a Student-t
Rt = ones(NbD,NSim)*m; % initialize at mean
Sigt2=zeros(NbD,NSim); % place holder for variance
Sigt2(1,:)=(sig0^2).*ones(1,NSim);

xt2=xt.^2;
for  i = 2:NbD
    ep = xt(i-1,:)> 0;
    en = xt(i-1,:)<=0;
    Sigt2(i,:) = w + ap.*ep.*xt2(i-1,:) + an.*en.*xt2(i-1,:) + Sigt2(i-1,:)*be; 
end

Rt = Rt + sqrt(Sigt2).*xt;
ST = S0.*exp(sum(Rt/100) ); % get time T and don't forget the scaling of the data
hist(ST,20)

dz=10;
z=3800:dz:5100; z=z';
kernz = kern(ST,-.65,z);
kernz = kernz/(sum(kernz)*dz);
plot(z,kernz)

subD=kernz;
save subD;