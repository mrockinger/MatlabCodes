function [beta,stderr,vc,logl]=Max_lik(lik_fct,b0,vc_type,A,b,Aeq,beq,lb,ub,nonlcon,options,varargin);
% Max_lik.m
% computes the maximum-likelihood estimates and associated standard errors
% allows for Sandwiched standard errors (White) or usual hessian ones
% implementation by Michael Rockinger March 20, 2004.
% please, report bugs to the author
% implemented own Hessian computation on 21.08.2005.
% it is possible that for large parameter values the Hessian behaves ugly

f0=feval(lik_fct,b0,varargin{:});
T=size(f0,1);

if T==1;
    error('Likelihood function should return a column vector of scores')
end

[beta,fval,exitflag,output,lambda,grad,hessian] =...
    fmincon(@ml_sum,b0,A,b,Aeq,beq,lb,ub,nonlcon,options,lik_fct,varargin{:});

hessian=HessMp(@ml_sum,beta,lik_fct,varargin{:});

inv_h=-hessian\eye(size(hessian,1));

if strcmp(vc_type,'Sandwich');
    
    disp('estimation of variance-covariance matrix via Sandwich [White]');
    g  = gradp(lik_fct,beta,varargin{:});
    vc = inv_h*(g'*g)*inv_h;
    
else %default is information matrix
    
    disp('estimation of variance-covariance matrix via Hessian'); 
    vc = -inv_h;
    
end

stderr  = sqrt(diag(vc));
logl    = -fval;
tstudent= beta./stderr;
parno   = (1:size(beta,1))';

if T-size(beta,1)>0;
     pvalue  = 2*(1-tcdf( abs(tstudent), T-size(beta,1) ));
 else
     pvalue = NaN*zeros(size(beta));
end

Res     = [ parno beta stderr tstudent pvalue];

fprintf('\n\n\n **********************************************************************\n');
if T-size(beta,1)<=0;
    fprintf('\nWarning\n')
    fprintf('Model contains more parameters than observations \n')
        fprintf('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n')
end
fprintf('Number of observations: %12.4f\n',T);
fprintf('Value of likelihood     %12.4f\n',-fval);
fprintf('Number of parameters    %12.4f\n',size(beta,1));
fprintf(' **********************************************************************\n');
fprintf('       parameter       beta        stderr    t-student      p-value\n');
fprintf('  %12.0f %12.4f  %12.4f %12.4f %12.4f\n', Res' );




function l=ml_sum(b,lik_fct,varargin)
l=feval(lik_fct,b,varargin{:});
l=-sum(l);

function g=gradp(f,x0,varargin)
% computes the gradient of f evaluated at x
% uses forward gradients. Adjusts for possible differently scaled x by taking percentage increments
% this function is the equivalent to the gradp function of Gauss
% f should return either a scalar or a column vector
% x0 should be a column vector of parameters
f0=feval(f,x0,varargin{:}); 
[T,c]=size(f0);

if size(x0,2)>size(x0,1)
    x0=x0';
end
k=size(x0,1); % number of parameters wrt which one should compute gradient

h=0.0000001; %some small number

g=zeros(T,k); %will contain the gradient
e=eye(k); 
for j=1:k;
    if x0(j)>1; % if argument is big enough, compute relative number   
        f1=feval(f,(x0.*( ones(k,1) +  e(:,j) *h )),varargin{:});    
        g(:,j)=(f1-f0)/(x0(j)*h);    
    else
        f1=feval(f, x0 +  e(:,j) *h ,varargin{:});    
        g(:,j)=(f1-f0)/h;    
    
    end
    
end

function H=HessMp(f,x0,varargin)
% computes the Hessian matrix of f evaluated at x0. If x0 has K elements,
% function returns a KxK matrix. The function f, given the ML context, is
% expected to be a column vector
% uses central differences 

% f should return either a scalar or a column vector
% x0 should be a column vector of parameters
f0=feval(f,x0,varargin{:}); 
[T,co]=size(f0);
if co>1; error('Error in HessMp, The function should be a column vector or a scalar'); end

[k,c]=size(x0);
if k<c,
    x0=x0';
end
k=size(x0,1); % number of parameters wrt which one should compute gradient

h=0.00001; %some small number

H=zeros(k,k); %will contain the Hessian
e=eye(k); 

h2=h/2;
for ii=1:k;
      if x0(ii)>100; % if argument is big enough, compute relative number   
        x0P= x0.*( ones(k,1) +  e(:,ii) *h2 );
        x0N= x0.*( ones(k,1) -  e(:,ii) *h2 );
        Deltaii = x0(ii)*h;
    else
        x0P = x0 +  e(:,ii) *h2;
        x0N = x0 -  e(:,ii) *h2;
        Deltaii = h;
    end
    
    for jj=1:ii
    if x0(jj)>100; % if argument is big enough, compute relative number   
        x0PP = x0P .* ( ones(k,1) +  e(:,jj) *h2 );
        x0PN = x0P .* ( ones(k,1) -  e(:,jj) *h2 );
        x0NP = x0N .* ( ones(k,1) +  e(:,jj) *h2 );
        x0NN = x0N .* ( ones(k,1) -  e(:,jj) *h2 );
        Delta = Deltaii*x0(jj)*h;
    else
        x0PP = x0P  +  e(:,jj) *h2; 
        x0PN = x0P  -  e(:,jj) *h2; 
        x0NP = x0N  +  e(:,jj) *h2; 
        x0NN = x0N  -  e(:,jj) *h2; 
        Delta = Deltaii*h;
    end
    
        fPP = feval(f,x0PP,varargin{:});   % forward,forward
        fPN = feval(f,x0PN,varargin{:});   % forward,backward
        fNP = feval(f,x0NP,varargin{:});    % backward,forward
        fNN = feval(f,x0NN,varargin{:});    % backward,backward
        
        H(ii,jj)=(sum(fPP)-sum(fPN)-sum(fNP)+sum(fNN))/Delta;
        H(jj,ii)=H(ii,jj);
    end
end