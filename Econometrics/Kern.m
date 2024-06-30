function f = kern(x,h,b);
% the kernel estimator
% this is a Matlab translation of Bruce Hansen's kernel estimator
npt = rows(b);
f   = zeros(npt,1); %will contain the kernel estimate

if h<=0
     if h==0
          c=1; 
     else
          c=-h
     end
            
     sx = std(x);
     h  = c*(npt^(-1/5))*sx;
end
     
n = rows(x);
for j=1:npt;
       f(j) = sum(pdfn((b(j)-x)/h))/(n*h);
end
     
function y=pdfn(x)
y=exp(-0.5*x.^2)/sqrt(2*pi);