function [y_simul,a_simul] = Kalman_Simulation(Z,d,T,c,R,a_star,H,Q,timevar,ns)
% simulation of a state-space model
% code is a Matlab translation of Thierry Roncalli's Gauss code.


  at=a_star; m=rows(at);
  if timevar == 1;
    n=cols(Z)/m; g=cols(R)/m;
  else;
    n=rows(Z); g=cols(R);
  endif;

  if timevar /= 1;
    Zt=Z;
    dt=d;
    Ht=H;
    Tt=T;
    ct=c;
    Rt=R;
    Qt=Q;
    if Ht == zeros(n,n)
      P1=0;
    else
      P1=chol(Ht)';
    end
  P2=chol(Qt)';
  endif;

  % Simulation 

  y_simul=zeros(ns,n);
  a_simul=zeros(ns,m);

  i=1;
  do until i>ns;

    if timevar == 1;
      Zt=reshape(Z(i,:),n,m);           % Z(t) 
      dt=d(i,:)';                       % d(t) 
      Ht=reshape(H(i,:),n,n);           % H(t) 
      Tt=reshape(T(i,:),m,m);           % T(t) 
      ct=c(i,:)';                       % c(t) 
      Rt=reshape(R(i,:),m,g);           % R(t) 
      Qt=reshape(Q(i,:),g,g);           % Q(t) 
      if Ht == zeros(n,n);
        P1=0;
      else;
        P1=chol(Ht)';
      endif;
      P2=chol(Qt)';
    endif;

    epsilon=P1*rndn(n,1);
    nu=P2*rndn(g,1);

    at = Tt*at + ct + Rt*nu;
    yt = Zt*at + dt + epsilon;

    a_simul(i,:) = at';
    y_simul(i,:) = yt';

    i=i+1;
  endo;

 
 
