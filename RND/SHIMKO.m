function SHIMKO()
% Shimko.m 
% procedure to compute RNDs following Shimko.
% code follows Sophie Coutant, a former PhD of mine (great girl!)
% all responsibility for mistakes declined by authors

safe={}; 
G_SIGMA=0; 
G_K=0; 
_APPROXIMATION_COEFFICIENTS=0;

%*******************************************************************
%                      ***** loads data *****
%*******************************************************************


%*******************************************************************
%  **** Gives number of maturities and extracts information *****
%*******************************************************************

  for i = 1:nbmat

    nbase=selif(basea,basea[.,4].==matjour[i]);

    F = nbase[1,3]; 
    T= nbase[1,4];  
    KCall = nbase[.,5]; 
    rFr = nbase[1,6]; 
    PCall = nbase[.,7]; 
    S0=nbase[1,2]; 
    r0=nbase[1,12]; 
    nbK=rows(KCall); 
    smile=nbase[.,11];
    
    
    div=rfr-ln(F./S0)./T; 
    TF=0; 
    Fwd=Fwd|F;

     _m=2;

     smileshimko=Impvol(KCall,_m);
     
    { vnam,m,b,stb,vc,stderr,sigma,cx,rsq,resid,dwstat } = 
                  OLS("",smile,KCall~(KCall^2));

    b0=0.5;
    {bf,ub,gb,retcodeb}=optmum(&optimisat,b0);
    bf=(atan(bf)+pi/2)/pi;
 
    /*Shimko RND*/
    N=2001; Kmin=1500; Kmax=5000; h=(Kmax-Kmin)/N; z=seqa(KMin,h,N);
    coeff=b;
    dens=shim_dens(coeff,z); 
    dens=dens.*(dens.>0);
    shimko=shimko~dens;

    /*computing moments*/
    Esp=xtrapeze(z,z.*dens);
    Var=xtrapeze(z,((z-Esp)^2).*dens);
    Skew=xtrapeze(z,((z-Esp)^3).*dens)/(Var^1.5);
    Kurt=xtrapeze(z,((z-Esp)^4).*dens)/(Var^2);
    moments=(esp~var~skew~kurt);

    /*computing MSE and ARE*/
    prix = (American_Call_Shimko(F,T,Kcall,bf,smileshimko,rfr));
    MSE=sumc(((PCall-prix)^2))/(nbK-1);
    v=((PCall-prix)^2)./(PCall^2);
    ARE=sumc(v)/(rows(KCall)-1);
  
    % store resultats:date,nbmat,forward,mat,nb_strikes,volume, 
    %mean,standard deviation,Skewness,kurtosis,b3,b4,STE */
    
safe=safe|(stof(dat)~i~F~(T*365)~rows(KCall)~Esp~sqrt(Var)~Skew~Kurt~rFr~ARE~MSE
~coeff[1]~coeff[2]*1e4~coeff[3]*1e7~bf);

    i=i+1;
  endo;

% make some nice graphs with smiles
plot(shimko)
title("RNDs");
xlabel("F]T["); 
ylabel("Densite");

% *******************************************************************************

function sm=Impvol(u,m)

sm=approximation_polynomiale(KCall,smile,u,ones(rows(KCall),1),m).*(u.>=minc(KCa
ll) .and u.<=maxc(KCall))
	+ smile[1].*(u.<minc(KCall)) + smile[rows(KCall)].*(u.>maxc(KCall));
    
% ===============================================================================

function dist = optimisat(b)
b=(atan(b)+pi/2)/pi;
prix = (American_Call_Shimko(F,T,Kcall,b,smileshimko,rfr));
dist=PCall-prix;
dist=dist'dist;

% ================================================================================

function C=Call_Black(F0,T,K,sig,r)
d0 = ( ln(F0./K) + (0.5*sig^2).*T )./(sig.*sqrt(T));
d1 = d0 - sig.*sqrt(T);
C  = exp(-r.*T).* ( F0.*cdfn(d0) - K.*cdfn(d1) );

% ==============================================================================

function prix=call_Shimko(F,T,K,r,m)
vol=impvol(K,m);
prix=call_black(F,T,K,vol,r);

% ===============================================================================

function dCall = American_call_Shimko(F,T,K,b,sig,r)
w=b;
erT=exp(-r*T);
dist=call_Black(F,T,K,sig,r)/ert;
Cu=dist;
Cl=maxc( ((F-K)')|(erT*dist') );
dCall=w*Cu+(1-w)*Cl;

% ================================================================================

function prix = call_Shimko2(F,T,K,r,b)
vol=b[1]+b[2]*K+b[3]*K^2;
prix=call_black(F,T,K,vol,r);

% ================================================================================

function dCall=American_call_Shimko2(F,T,K,b,r)
w=b[4];
erT=exp(-r*T);
dist=call_Shimko2(F,T,K,r,b[1:3])/ert;
Cu=dist;
Cl=maxc( ((F-K)')|(erT*dist') );
dCall=w*Cu+(1-w)*Cl;

% ===============================================================================

function density = shim_dens(b,x)
st=sqrt(T);
a0=b[1];  a1=b[2];  a2=b[3];
sig=(a0+a1*x+a2*x^2).*(x.>=minc(KCall) .and x.<=maxc(KCall)) + 
smile[1].*(x.<minc(KCall)) + smile[rows(KCall)].*(x.>maxc(KCall));
v=sig*st;
ds=(a1+2*a2*x).*(x.>=minc(KCall) .and x.<=maxc(KCall));
dss=(2*a2).*(x>=minc(KCall) .and x<=maxc(KCall));

d1=(ln(F/x)+0.5*v^2)./v;
d2=d1-v;

d1x=-1/(sig.*st.*x) + 0.5*ds.*st + ln(F/x).*(-ds./(sig^2.*st));
d2x=d1x-ds.*st;

%d1xx=-(2*a2*st*v-2*(a1+2*a2*x)^2*T)./(v^3).*ln(F/x);
d1xx=d1xx + (a1+2*a2*x)*st./(x.*v^2);
d1xx=d1xx + (v+x*st.*(a1+2*a2*x))./(x.*v)^2+a2*st;@
d1xx=ds./((sig^2).*st.*x) + 1/(sig.*st.*(x^2));
d1xx=d1xx + 0.5*dss.*st + ds./((sig^2).*st.*x);
d1xx=d1xx - ln(F/x).*dss./((sig^2).*st);
d1xx=d1xx + 2*ln(F/x).*(ds^2)./(sig^3)/st;
d2xx=d1xx-dss*st;

n1=pdfn(d1);
n2=pdfn(d2);

density=(F.*n1.*(d1xx-d1.*d1x^2));
density=density - (n2.*(2*d2x+x.*d2xx-x.*d2.*d2x^2));

% ========================================================================

function i=xtrapeze(x,y)
n=rows(x);
i=0.5*sumc ( (y[2:n,.]+y[1:n-1,.]) .* (x[2:n,.]-x[1:n-1,.]) );

% ========================================================================

function t=approximation_polynomiale(x,y,z,poids,n)
  a=zeros(n+1,n+1);
  b=zeros(n+1,1);
  j=1;
  do until j>n+1;
    i=seqa(1,1,n+1);
    h=2*(n+1)-(i'+j);
    q=poids.*(x^h);
    a[j,.]=sumc(q)';
    b[j]=sumc(poids.*y.*(x^(n+1-j)));
    j=j+1;
  endo;
  coefficients=b/a;
  p=rows(z);
  zz=ones(p,1)~polymat(z,n);
  t=zz*rev(coefficients);
  _approximation_coefficients = rev(coefficients);

% ========================================================================

function y=rect(x,a,b);
ex=exp(x);
y=a + (ex./(1+ex)) * (b-a);


% ========================================================================

proc (1) = Call_Black(F0,T,K,sig,r);
local d0,d1,C;
  d0 = ( ln(F0./K) + (0.5*sig^2).*T )./(sig.*sqrt(T));
  d1 = d0 - sig.*sqrt(T);
  C = exp(-r.*T).* ( F0.*cdfn(d0) - K.*cdfn(d1) );
retp(C);
endp;


proc (1) = Put_Black(F0,T,K,sig,r);
local d0,d1,C;
  d0 = ( ln(F0./K) + (0.5*sig^2).*T )./(sig.*sqrt(T));
  d1 = d0 - sig.*sqrt(T);
  C = exp(-r.*T).*( -F0.*cdfn(-d0) + K.*cdfn(-d1) );
retp(C);
endp;

