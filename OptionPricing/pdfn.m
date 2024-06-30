function p=pdfn(d);
%pdfn returns the expression of the density of a standard normal

p = 1/sqrt(2*pi)*exp(-0.5*d.^2);