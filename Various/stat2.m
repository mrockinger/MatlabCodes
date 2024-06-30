% Stat2.M
% Simulates draws from a normal density and an exponential and 
% verify that the mean satisfies the unbiasedness property
%

NSim=100; % Number of samples to be drawn
MSim=10;  % size of each sample
lda=0.5;
mu=2;

mean_norm=zeros(NSim,1);
mean_expo=zeros(NSim,1);

for j=1:NSim

    % draw from exponential
    u=rand(NSim,1);
    expo=-1/lda*log(u);

    % draw from normal
    norm=mu+randn(NSim,1);

	mean_norm(j)=mean(norm);	
	mean_expo(j)=mean(expo);	

end

s=1:NSim; s=s';

mean_norm=cumsum(mean_norm)./s;
mean_expo=cumsum(mean_expo)./s;

plot(s,mean_norm,'-v',s,mean_expo,'-h');
title('Unbiasedness of mean');

