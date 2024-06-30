function TestSimPoiss()
% tests the Poisson generator
x=SimPoiss(3,1000) % 1000 draws of a poisson with mean 3
hist(x); 
fprintf('The mean should be equal to 3 %5.2f \n',mean(x));
