function cdfn=cdfn(x)
% cdfn cumulative normal density function
% builds on error function
% this is not a vectorized function

cdfn=0.5+0.5*erf(x/sqrt(2));
