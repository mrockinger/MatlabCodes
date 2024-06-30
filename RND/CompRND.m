function CompRNDFTSE()
% loads and compares various risk neutral densities
clc; cla;
z=3800:10:5100; z=z'; %support for RND
load BenchRND;
load MixRND;
load GB2RND;
load SempRND;
%load ShimkoRND;
load subD; % subjective density at a 20 days time horizon estimated via a kernel fit

for i=1:1;
   hold on
   plot(z,BenchRND(:,i),'k','LineWidth',2);
   plot(z,MixRND(:,i),'g');
   plot(z,GB2RND(:,i),'r');
 %  plot(z,ShimkoRND(:,i),'c');
   plot(z,SempRND(:,i),'k');
 
   plot(z,subD,'m');
   hold off
end

