function TempAgreg();
% investigates in terms of moments the effects of temporal aggregation
clc;

load rtSP;
Rt=rtSP;
size(Rt)

%plot(Rt)
%title('SP500 fréquence journalière')
%[mi,vmi]=min(Rt);

%Rt=Rt(4640:4660)

momM = zeros(240,4);

%rt=rt(1:23)

for i=1:240;
    
    x=getx(Rt,i); % agregate data up to a given frequency
    
   % plot(abs(x));
   % title(' Val. Abs. SP500 fréquence journalière');
   % return
    m=getmom(x);
    
    %scale the variables
    m(1) = m(1)*(250/i);
    m(2) = m(2)*sqrt(250/i);
    
    momM(i,:)=m;
end

figure
plot(momM(:,1));
title('Moyenne annualisée pour des fréquences grandissantes');

figure
plot(momM(:,2));
title('Volatilité annualisée pour des fréquences grandissantes');

figure
plot(momM(:,3));
title('Skewness pour des fréquences grandissantes');

figure
plot(momM(:,4));
title('Kurtosis pour des fréquences grandissantes');



function x=getx(rt,i)
% grabs data for a given frequency.
T  = rows(rt);
nT = floor(T/i);
x  = reshape(rt(1:nT*i),i,nT);
x  = sum(x,1)';

function m=getmom(x)
% here compute the various moments
m    = zeros(1,4);
m(1) = mean(x);
m(2) = std(x);
m(3) = skewness(x);
m(4) = kurtosis(x);