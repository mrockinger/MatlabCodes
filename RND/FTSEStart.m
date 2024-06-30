function FTSEStart()
% in this file load option data. Some code that should get you going.
% the chosen date is March 26, 2004.
clc;

% load data
A=xlsread('FTSE100Mar2604.xls');
%A

%Reshape the data
Ft  = A(2,1); % FT100
Kv  = A(1,2:17); % strikes
Tv  = A(2:6,18)/365; % days till maturity
rtv = log(1+A(2:6,19)/100); % the interest rate in continous time
rtv
%return

CP  = A(2:6,2:17); % all the call and put prices

M   = size(CP,1); % number of dates
N   = size(CP,2); % number of options for one date (both C+P)

Cidx=1:2:N;
Pidx=2:2:N;
Cv=zeros(M*N/2,1);
Pv=zeros(M*N/2,1);
for i=1:5;
    for j=1:N/2;
        Cv(1+(i-1)*N/2+(j-1))=CP(i,1+2*(j-1));
        Pv(1+(i-1)*N/2+(j-1))=CP(i,2+2*(j-1));
    end
end

Opt = [Cv; Pv];

S0  = Ft(1,1)*ones(size(Opt,1),1);
K   = Kv(2:2:N)';

K   = kron(ones(5,1),K);
r   = kron(rtv,ones(size(Cv,1)/5,1));
T   = kron(Tv,ones(size(Cv,1)/5,1));

CPind = [ones(size(Cv,1),1) ;zeros(size(Cv,1),1)];
AllInfo= [S0 [K; K] CPind [r; r] [T;T] [Cv; Pv]];

% control of data
for i=1:size(AllInfo,1);
   fprintf('%8.2f %8.2f %3.0f %8.2f %8.2f %8.2f\n', AllInfo(i,:));
end

%
% From now on you are on your own. Good luck!
%