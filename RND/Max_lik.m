function [beta,stderr,vc,logl]=Max_lik(lik_fct,b0,vc_type,A,b,Aeq,beq,lb,ub,nonlcon,options,varargin);
% Max_lik.m
% computes the maximum-likelihood estimates and associated standard errors
% allows for Sandwiched standard errors (White) or usual hessian ones
% implementation by Michael Rockinger March 20, 2004.
% please, report bugs to the author

f0=feval(lik_fct,b0,varargin{:});
T=size(f0,1);

[beta,fval,exitflag,output,lambda,grad,hessian] =...
    fmincon(@ml_sum,b0,A,b,Aeq,beq,lb,ub,nonlcon,options,lik_fct,varargin{:});

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
pvalue  = 2*(1-tcdf( abs(tstudent), T-size(beta,1) ));
Res     = [ parno beta stderr tstudent pvalue];


fprintf('\n\n\n **********************************************************************\n');
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
        
    f1=feval(f,(x0.*( ones(k,1) +  e(:,j) *h )),varargin{:});    
    g(:,j)=(f1-f0)/(x0(j)*h);    
    
end
