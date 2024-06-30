function g = GTdens(z,lda,eta);
% constructs Hansen's generalized Student t
ep1 = (eta+1)/2;
c   = gamma(ep1)/(sqrt( pi * (eta-2) )*gamma( eta/2 ));
a   = 4*lda*c*(eta-2)/(eta-1);
b   = sqrt(1 + 3*lda^2 - a^2);
arg = b*z+a;
ind = arg<0;
gn  = b.*c.*(1 + 1./(eta-2) * ( arg./(1-lda) ).^2 ).^(-ep1); %negative part of density
gp  = b.*c.*(1 + 1./(eta-2) * ( arg./(1+lda) ).^2 ).^(-ep1);
g   = gn.*ind + gp.*(1-ind);

