function [mp,sp]=min_var_front(m1,s1,m2,s2,rho, delta);
% returns the expected return-risk caracteristics
% for a portfolio containint delta% of portfolio 1 and 1-delta
% of portfolio 2
mp=delta.*m1+(1-delta).*m2;
omd=-(delta-1);
sp1=(delta.^2).*(s1.^2)+(omd.^2).*(s2.^2)+delta.*omd.*s1.*s2.*rho.*2;
sp=sqrt(sp1);
