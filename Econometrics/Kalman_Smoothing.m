function [a_smooth,P_smooth] = Kalman_Smoothing(y,Z,d,T,c,R,a0,P0,H,Q,timevar);
% Kalman smoothing
% code is a Matlab translation of Thierry Roncalli's.

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

a_cond   = zeros(nobs,m); 
a        = zeros(nobs,m);
P_cond   = zeros(nobs,m*m); 
P        = zeros(nobs,m*m);
a_smooth = zeros(nobs,m); 
P_smooth = zeros(nobs,m*m);

if timevar ~= 1 % then matrices are time invariant
  Zt=Z;
  dt=d;
  Ht=H;
  Tt=T;
  ct=c;
  Rt=R;
  Qt=Q;
end

for  i=1:nobs

    yt=y(i,:)';                         % y(t) 

    if timevar == 1
      Zt=reshape(Z(i,:),n,m);           % Z(t) 
      dt=d(i,:)';                       % d(t) 
      Ht=reshape(H(i,:),n,n);           % H(t) 
      Tt=reshape(T(i,:),m,m);           % T(t) 
      ct=c(i,:)';                       % c(t) 
      Rt=reshape(R(i,:),m,g);           % R(t) 
      Qt=reshape(Q(i,:),g,g);           % Q(t) 
    end;

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

    % Save results 

    a(i,:)=at';
    a_cond(i,:)=at_1';
    P(i,:)=vecr(Pt)';
    P_cond(i,:)=vecr(Pt_1)';

end;

% Smoothing Equations 

a_smooth(nobs,:)=at';
P_smooth(nobs,:)=vecr(Pt)';

for i=(nobs-1):-1:1;
    
    if timevar == 1;
      Tt=reshape(T(i+1,:),m,m);           % T(t+1) 
    end;

    Pt=reshape(P(i,:),m,m);               % P(t)     
    Pt_1=reshape(P_cond(i+1,:),m,m);      % P(t+1|t) 
    Pt_1_T=reshape(P_smooth(i+1,:),m,m);  % P(t+1|T) 
    at=a(i,:)';                           % a(t)     
    at_1=a_cond(i+1,:)';                  % a(t+1|t) 
    at_1_T=a_smooth(i+1,:)';              % a(t+1|T) 

    inv_Pt_1 = Pt_1\eye(size(Pt_1,1));         % Inversion de Ft           

    P_star = Pt*Tt'*inv_Pt_1;

    as = at + P_star*(at_1_T-at_1) ;          % a(t|T) formula(3.6.16a) 
    Ps = Pt + P_star*(Pt_1_T-Pt_1)*P_star' ;  % P(t|T) formula(3.6.16b) 

    a_smooth(i,:)=as';
    P_smooth(i,:)=vecr(Ps)';
    
end


