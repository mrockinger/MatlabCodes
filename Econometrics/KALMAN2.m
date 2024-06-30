function Kalman2()
% some more Kalman filtering testing program
% time varying parameter model
% code is a transcription of Thierry Roncalli's Gauss code.
nobs=200;
b0=recserar(randn(nobs,1)*1.5,10,1);
b1=recserar(randn(nobs,1)*0.2,4,1);

%plot([b0 b1]);
%return

x0 = ones(nobs,1); 
x1 = randn(nobs,1)*25; 
u  = randn(nobs,1)*2;
y  = x0.*b0 + x1.*b1 + u;
X  = [x0 x1]; 
beta_mco=(X'*X)\(X'*y)
%return
%plot([y x1])
%return

Z=[x0 x1]; 
d=zeros(nobs,1); 
c=zeros(nobs,2); 
T=[ 1 0; 0 1 ];
T=kron(vecr(T)',ones(nobs,1)); 
R=T;

a0=0;
P0=0;
H=0;
Q=0;

beta0=ones(5,1);
lb=[-10 -10 0.00001 0.0001 0.0001]; % sigma_epsilon sigma_eta
ub=[+10 +10 100 100 100];
lb=lb';
ub=ub';

options=optimset('Diagnostics','on','Display','iter');
[beta,stderr1,vc,logl]=Max_lik(@ml,beta0,'Hessian',[],[],[],[],lb,ub,[],options,...
           y,Z,d,T,c,R,a0,P0,H,Q,1);

 beta(3:5)=beta(3:5).^2;
 a0= beta(1:2); 
 P0=[ beta(4) 0 ; 0 beta(5) ];
 H = kron(beta(3),ones(nobs,1)); 
 Q= kron(vecr(P0)',ones(nobs,1));
 
 [y_cond,v,a,a_cond,P,P_cond,F,logl]=kalman_filter(y,Z,d,T,c,R,a0,P0,H,Q,1);
 [a_smooth,P_smooth]=kalman_smoothing(y,Z,d,T,c,R,a0,P0,H,Q,1);

s=1:nobs;
subplot(2,1,1);
plot(s,a_cond,s,a_smooth,s,b0,s,b1);

subplot(2,1,2);
plot(s,P_cond,s,P_smooth);

disp('Hit a key to continue');
pause;


%forecast over half of sample
[y_prev,a_prev,P_prev,MSE] = Kalman_Forecasting(y,Z,d,T,c,R,a0,P0,H,Q,1,floor(size(y,1)/2));

subplot(2,1,1);
plot(s,y,s,y_prev)

subplot(2,1,2);
plot(s,a,s,a_prev)


%*****************************************************

function logl=ml(beta,y,Z,d,T,c,R,a0,P0,H,Q,xx);
  beta(3:5)=beta(3:5).^2;
  a0=beta(1:2); 
  P0=[ beta(4) 0; 0 beta(5) ];
  nobs=size(y,1);
  H=kron(beta(3),ones(nobs,1)); 
  Q=kron(vecr(P0)',ones(nobs,1));
  [y_cond, v, a, a_cond, P, P_cond,F,logl]=kalman_filter(y,Z,d,T,c,R,a0,P0,H,Q,xx);
