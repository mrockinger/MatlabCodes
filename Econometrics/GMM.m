function [ beta, stderr, covbeta, Qmin, test, ptest ] = GMM(MomFct,beta0,...
                A,b,Aeq,beq,lb,ub,nonlcon,options,...
                GMMLags,GMMiter,GMMtol1,GMMtol2,varargin);
% an interface for GMM estimation
% also prints the final parameter with standard errors
% The GMM code comes from T. Roncalli.

% start with a call to the core function to get an idea of the sample size involved
fs=feval(MomFct,beta0,varargin{:});
[T,r]=size(fs); %T=number of observations

k=size(beta0,1); %number of parameters
ddl=T-k; %degree of freedom

b0=beta0;
Omega=eye(r); %start with identity matrix 
Qmino=100000; %something big
i=1;
while i <= GMMiter;
    
      betai=b0; % parameters we started with 
           
      [beta,Qmin,exitflag,output,lambda,grad,hessian] =...
      fmincon(@GMM_obj,b0,A,b,Aeq,beq,lb,ub,nonlcon,options,Omega,MomFct,varargin{:});
      disp('parameters after optimization');
      beta
    
      if ((max(abs(beta-b0))<GMMtol1) | (abs(Qmin-Qmino)<GMMtol2)) %complete with precision over J
       %   fprintf('the parameters before and after  loop are');
       %   [betai beta ]
          break;
      end
      
      H  = feval(MomFct,beta,varargin{:});
      mH = mean(H);
      H  = H - kron(mH,ones(size(H,1),1));
      Omega = (H'*H)/T;
   
      for j=1:GMMLags
        % size(H(j+1:end,:)')
        % size(H(1:end-j))
          Gamma = ( H(j+1:end,:)'*H(1:end-j,:)) /T;
          Omega = Omega + (1-j/(GMMLags+1)) * (Gamma + Gamma'); %Parzen type weighting
      end
   
      b0=beta;
      fprintf('the weighting matrix is \n');    
      Omega
      fprintf('done with GMM iteration %5.0f\n',i);
      Qmino=Qmin;
    
      i=i+1;
end

D=gradp(@GMM_momgen,beta,MomFct,varargin{:});

if cond(Omega)>100000
    invOmega=pinv(Omega); %Moore-Penrose inverse
else
    invOmega=Omega\eye(size(Omega,1));
end

DiOmegaD=D'*invOmega*D;
covbeta=((DiOmegaD)\eye(size(DiOmegaD,1)))/T;
stderr=sqrt(diag(covbeta));
corbeta=covbeta./kron(stderr,stderr');
tstudent=beta./stderr;
pvalue=2*(1-tcdf(abs(tstudent),ddl));

r=size(Omega,1);
k=size(DiOmegaD,1);
if r>k
    test=T*Qmin;
    ptest=chi2cdf(test,r-k);
else
    test=[];
    ptest=[];
end

name=1:size(beta,1);
name=name';
omat=[ name beta stderr tstudent pvalue];
fprintf('\n');
fprintf('Number of observations : .....................%12.4f\n',T);
fprintf('Number of parameters :   .....................%12.4f\n',k);
fprintf('Number of degrees of freedom : ...............%12.4f\n',ddl);
fprintf('Number of orthogonalit conditions :...........%12.4f\n',r);
fprintf('Value of the objective function : ............%12.4f\n', Qmin);
fprintf('\nTest of overidentification of restrictions :..%12.4f \n',test);
fprintf('\nCorresponding marginal probability : .........%12.4f \n',ptest);
fprintf('\n');
fprintf('\n');
fprintf(' Parameter  Estimate      Standard Error   Student-t    Signif.\n');
fprintf('--------------------------------------------------------------------\n');
for j=1:size(omat,1);
      fprintf('%7.0f %14.6f %14.6f %14.6f %14.6f \n',omat(j,:));
end
fprintf('\n');

fprintf('correlation matrix of parameters\n');
for j=1:size(corbeta,1);
      fmt1='%8.4f';
      fmt=kron(ones(1,size(corbeta,1)),fmt1);
      fprintf([fmt '\n'],corbeta(j,:));
end

%----------------------------------------------------------------

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

h=0.000000000001; %some small number

g=zeros(T,k); %will containt the gradient
e=eye(k); 
for j=1:k;
        
    f1=feval(f,(x0.*( ones(k,1) +  e(:,j) *h )),varargin{:});    
    g(:,j)=(f1-f0)/(x0(j)*h);    
    
end

%----------------------------------------------------------------

function GT=GMM_momgen(beta,fnam,varargin);
% computes the average over momentized observations. Corresponds
% to computing gT.
% imports the name of the function, the parameters and the observations

mt=feval(fnam,beta,varargin{:});
GT=mean(mt)';

%----------------------------------------------------------------

function J=GMM_obj(beta,Omega,fnam,varargin);
% computes the value of the objective function.
% imports the name of the function, the parameters and the observations

GT = feval(@GMM_momgen,beta,fnam,varargin{:});
if cond(Omega)>100000
    invOmega=pinv(Omega); %Moore-Penrose inverse
else
    invOmega=Omega\eye(size(Omega,1));
end
J=GT'*invOmega*GT;

%----------------------------------------------------------------
