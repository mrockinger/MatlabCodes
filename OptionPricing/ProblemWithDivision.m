function ProblemWithDivision()
% shows how studpid Matlab is with divisions. Here the programmer forgot to
% type ./ and typed just /.
St=100*exp(cumsum(1+2*randn(10,1)));
rt=log(St(2:end)/St(1:end-1));
figure(1)
plot(rt)

figure(2)
rt=log(St(2:end)./St(1:end-1));
plot(rt)


