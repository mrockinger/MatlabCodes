% BinomCP.m
% tests the BinomTree.m function
% 
% Input: S0, K, T, sig, rho, N, CallPut
% CallPut>0 means compute call option. Else, compute Put option.
% Presumes European options
%
% Output: The tree with prices, the tree with option prices, the call or put option value
%
% Remark: T must be provided in yearly format. I.e. 60 days means T=60/365
%
%
clc;
S0=100;
K=100;
rho=log(1+0.08);
N=10;
sig=20/100;
T=60/365;
CallPut=1; %Positive corresponds to call and negative to get put price

BS=BSCall(S0,K,sig,T,rho);
fprintf('The price of Black-Scholes call is %9.4f\n',BS);

[StM,CtM,CP]=BinomTree(S0,K,rho,sig,T,N,CallPut);


fid=fopen('d:/finbox/option/out.out','w'); % for an output file
fid=1; % identifier for screen

fmt1='%12.4f';
fmt =Kron(ones(1,size(StM,2)),fmt1);
fmt=[fmt '\n'];
    
fprintf(fid,fmt,StM');
fprintf(fid,' \n\n');
fprintf(fid,fmt,CtM');

fprintf(fid,'\n\nthe price of the call option is %8.4f\n',CP);
fprintf(1,'\n\nthe price of the call option is %8.4f\n',CP);

if fid~=1;
    fid=fclose(fid);
end