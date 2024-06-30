function getFTSE()
% getFTSE.m
% extracts out of matrix FTSE data
% constructs returns and saves them
A=xlsread('StkIdx.xls');
FTSE=A(:,4);
FTSE(1:10)
FTSE(end-3:end)
rtFT=100*log( FTSE(2:end)./FTSE(1:end-1) );
rtFT(1:10,:)

subplot(2,1,1);
plot(FTSE);
title('FTSE 100 in level');
subplot(2,1,2);
plot(rtFT);
title('FTSE 100 returns');

save rtFT
