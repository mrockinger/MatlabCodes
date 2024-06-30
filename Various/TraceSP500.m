% Traces actual data


A = xlsread('StkIdx.xls');
% 1=date sp nik ftse dax cac
x=A(:,2);


h=plot(x);
set(h,'LineWidth',1);
title('SP500 1970-2001');
ylabel('SP 500 level');
xlabel('day');
%axis([0 1000 0 2500]);

