function [y_prev,a_prev,P_prev,MSE] = ...
         Kalman_Forecasting(y,Z,d,T,c,R,a0,P0,H,Q,timevar,t0);
% kalman_forecasting.m
% performs forecasting of Kalman filter
% code is a Matlab translation of Thierry Roncalli's Kalman filter code.
% forecast starts at t0

nobs = size(y,1); 
n    = size(y,2); 
at   = a0; 
at_1 = a0; 
Pt   = P0; 
Pt_1 = P0; 
logl = zeros(nobs,1);

if timevar == 1
    m=size(Z,2)/n; 
    g=size(R,2)/m;
else
    m=size(Z,2); 
    g=size(R,2);
end

y_prev=zeros(nobs,n)*NaN; 
a_prev=zeros(nobs,m)*NaN;
P_prev=zeros(nobs,m*m)*NaN; 
MSE   =zeros(nobs,n*n)*NaN;

if timevar ~= 1
    Zt=Z;
    dt=d;
    Ht=H;
    Tt=T;
    ct=c;
    Rt=R;
    Qt=Q;
end

for i=1:(t0-1)

    yt=y(i,:)';                         % y(t) 

    if timevar == 1
      Zt=reshape(Z(i,:),n,m);           % Z(t) 
      dt=d(i,:)';                       % d(t) 
      Ht=reshape(H(i,:),n,n);           % H(t) 
      Tt=reshape(T(i,:),m,m);           % T(t) 
      ct=c(i,:)';                       % c(t) 
      Rt=reshape(R(i,:),m,g);           % R(t) 
      Qt=reshape(Q(i,:),g,g);           % Q(t) 
    end

    % Prediction Equations 

    at_1 = Tt*at + ct ;                  % a(t|t-1) formula(3.2.2a) 
    Pt_1 = Tt*Pt*Tt' + Rt*Qt*Rt' ;       % P(t|t-1) formula(3.2.2b) 

    % Innovations 

    yt_1 = Zt*at_1 + dt ;                % y(t|t-1) formula(3.2.18) 
    vt = yt-yt_1 ;                       % v(t)     formula(3.2.19) 

    % Updating Equations 

    Ft = Zt*Pt_1*Zt' + Ht ;              % F(t)     formula(3.2.3c) 
    inv_Ft = Ft\eye(size(Ft,1));         % Inversion de Ft           

    at = at_1 + Pt_1*Zt'*inv_Ft*vt ;      % a(t)     formula(3.2.3a)   
    Pt = Pt_1 - Pt_1*Zt'*inv_Ft*Zt*Pt_1 ; % P(t)     formula(3.2.3b)   

end

  % Forecasting Equations 

for i=t0:nobs

    if timevar == 1
      Zt=reshape(Z(i,:),n,m);           % Z(t) 
      dt=d(i,:)';                       % d(t) 
      Ht=reshape(H(i,:),n,n);           % H(t) 
      Tt=reshape(T(i,:),m,m);           % T(t) 
      ct=c(i,:)';                       % c(t) 
      Rt=reshape(R(i,:),m,g);           % R(t) 
      Qt=reshape(Q(i,:),g,g);           % Q(t) 
    end

    at = Tt*at + ct ;                   % a(T+l|T) formula(3.5.7a)      
    Pt = Tt*Pt*Tt' + Rt*Qt*Rt' ;        % P(T+l|T) formula(3.5.7b)      
    yt = Zt*at + dt ;                   % y(T+l|T) formula(3.5.6a)      
    MSEt = Zt*Pt*Zt' + Ht ;             % MSE(y(T+l|T)) formula(3.5.6b) 

    a_prev(i,:) = at';
    y_prev(i,:) = yt';
    P_prev(i,:) = vecr(Pt)';
    MSE(i,:) = vecr(MSEt)';
end
