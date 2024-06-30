% Stat3.m
% Simulates draws from a normal density and an exponential and 
% checks what happens with the variance


NSim=1000; % Number of samples to be drawn
MSim=3;   % size of each sample
lda=1;
mu=2;

std_norm=zeros(NSim,1);
std_expo=zeros(NSim,1);

for j=1:NSim;

    % draw from exponential
    u=rand(MSim,1);
    expo=-1/lda*log(u);

    % draw from normal
    norm=mu+randn(MSim,1);

    % biased estimates
	std_norm(j)=sum( (norm-mean(norm)).^2 )/MSim;	
	std_expo(j)=sum( (expo-mean(expo)).^2 )/MSim;	

end;

s=1:NSim;s=s';

std_norm=cumsum(std_norm)./s;
std_expo=cumsum(std_expo)./s;

theovar=ones(NSim,1);
plot(s,std_norm,'-',s,std_expo,'--',s,theovar,'-');
title('Unbiasedness of standard error');
% axis([0 1 0.1 1]);