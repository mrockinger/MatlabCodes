function testkern()
% KERN.m
% Bruce Hansen's kernel code
% transcription from Gauss into Matlab 
% kern(x,h,b)
% 
% The pdf of the vector x is estimated non-parametrically by the kernel
% method, using pdf-N(0,1) as the kernel.
% A window of width h is used.
%   If h=0, default=npt^(-1/5)*std(x)         }
%   If h<0, default=npt^(-1/5)*std(x)*abs(h)  } see Singh and Ullah (1986)
%   If h>0, accept value without modification }
% 
% The points at which estimates are required are in vector b, which is nptx1.
% Vector b must be sorted.
%

% To test the function
z=randn(10000,1);
miz=min(z);
maz=max(z);
dz=((maz-miz)/1000);
zgrid=miz:dz:maz; zgrid=zgrid';
z=sort(z,1);
kernz = kern(z,-.65,zgrid);
kernz = kernz/(sum(kernz)*dz);

plot(zgrid,kernz);


