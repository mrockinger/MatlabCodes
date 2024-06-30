% illustrates graphical capabilities
x=randn(10,1);
plot(x,'--+'); % plots x against time
axis([0 13 -3 3]);
title('a nice simulation')
xlabel('time goes by')
ylabel('the realizations')
legend('name of curve')
text(8,-2,'Right lower corner');

% the following sets the color of all curves to black and changes style of
% lines to solid, dash-dot, dash-dash, and dotted
%set(0,'DefaultAxesColorOrder',[0 0 0],...
%      'DefaultAxesLineStyleOrder','-|-.|--|:')