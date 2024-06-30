function RRA()
% loads and compares various risk neutral densities
% computes the relative risk aversion
clc; cla;
dz =10;
z=3800:dz:5100; z=z'; %support for RND
load BenchRND;
load MixRND;

load GB2RND;
load subD; % subjective density at a 20 days time horizon estimated via a kernel fit
%load SempRND;
%load ShimkoRND;

for i=1:1;
   hold on
   p=subD; % objective density

   q=MixRND(:,i); % subjective density 
   [RRA1, RA1, zmi]=get_RA(z,p,q);
   
   q=GB2RND(:,i);
   [RRA2, RA1, zmi]=get_RA(z,p,q);
   
   plot(zmi,RRA1,zmi,RRA2);
      
   hold off
end

%-----------------------------------------------------------------------

function [RRA, RA, zmi]=get_RA(z,od,rd)
% constructs the coefficient of risk aversion and of relative risk aversion
logd = log(od./rd);                 % ratio of objective to risk-neutral density
RA   = diff(logd)/(z(2)-z(1));      % derivative equals risk aversion
zmi= ( z(1:end-1) + z(2:end) ) / 2; % this is midpoint of each z increment
RRA=zmi.*RA;                        % relative risk aversion
