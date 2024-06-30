% Simulates a sample with 1000 draws from a uniform (0,1)
% distribution
T=1000;
x=rand(T,1);

plot(x,'o');
set(findobj('Type','line'),'Color','k')
title('1000 draws from an uniform U(0,1)');
xlabel('index');
ylabel('u_n');