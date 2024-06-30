function TestUnifUnique()
% 
clc;
%UnifUnique(10,1,5);
disp('random shuffle of numbers')
UnifUnique(10,0,9)
disp('set of 10 unique indices out of a set of 20 possible ones')
UnifUnique(10,1,20)