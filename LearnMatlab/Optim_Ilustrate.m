function Optim_Ilustrate()
% illustration how to use the optimization package
clc;



options=optimset('MaxFunEvals', 1000, ...
'TolFun',1e-5, ...
'TolX',1e-5, ...
'Diagnostics','on', ...
'Display','iter');

beta0=[1 2]; %starting valuse
A=3;
[beta,Qmin,exitflag]=fminunc(@obj_func,beta0,options,A);
disp('the optimum is');
beta

disp('wheras it should be ');
[A;-1]

function z=obj_func(beta,A)
x=beta(1);
y=beta(2);
z=(x-A)^2 + (y+1)^2;

