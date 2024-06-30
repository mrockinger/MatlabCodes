function KLocalAR();
% estimates for the PURSE data set of Harvey an AR(1) with noise model
% code is a transcription of Thierry Roncalli's Gauss code.
clc;
echo on;

purse =[
   10    15    10    10    12    10     7    17    10    14     8    17    14 ...
   18     3     9    11    10     6    12    14    10    25    29    33    33 ...
   12    19    16    19    19    12    34    15    36    29    26    21    17 ...
   19    13    20    24    12     6    14     6    12     9    11    17    12 ...
    8    14    14    12     5     8    10     3    16     8     8     7    12 ...
    6    10     8    10     5     7 ];

% initialize data
y=purse'; 
Z=1; 
d=0; 
c=0; 
R=1; 
T=0; % just any value, nobody cares
a0=y(1); 
P0=0; 
H=0; 
Q=0;

lb=[0.0001; 0.0001; 0.5]; % sigma_epsilon sigma_eta auto-regressive parameter
ub=[100; 100; 1];
beta0=[std(y)^2; 1; 0.9]

options=optimset('Diagnostics','on','Display','iter');

[beta,stderr1,vc,logl]=Max_lik(@ml,beta0,'Sandwich',[],[],[],[],lb,ub,[],options,...
 y,Z,d,T,c,R,a0,P0,H,Q,0);

H=beta(1); Q=beta(2); T=beta(3);
[y_cond,v,a,a_cond,P,P_cond,F,logl] = kalman_filter(y,Z,d,T,c,R,a0,P0,H,Q,0);
w=v./sqrt(F);

subplot(2,1,1);
s=1:size(y,1);
plot(s,y,s,y_cond);

subplot(2,1,2);
plot(s,w);
title('Innovations standardisees');
disp('Hit a key to continue');
pause

%forecast over half of sample
[y_prev,a_prev,P_prev,MSE] = Kalman_Forecasting(y,Z,d,T,c,R,a0,P0,...
                   H,Q,0,floor(size(y,1)/2));

subplot(2,1,1);
plot(s,y,s,y_prev)

subplot(2,1,2);
plot(s,a,s,a_prev)

%********************************************

function logl=ml(beta,y,Z,d,T,c,R,a0,P0,H,Q,xx)
% evaluates the log-likelihood of model
H=beta(1); Q=beta(2); T=beta(3);
[y_cond,v,a,a_cond,P,P_cond,F,logl]=kalman_filter(y,Z,d,T,c,R,a0,P0,H,Q,xx);
