function showeval()
% illustrates how 'feval' operates
% here compute gradients for two trivial functions at trivial points

x=1.0;
c=2;
d=10;
grad(@f,x,c,d) % evaluates grandient of f at x

function y=f(x,d,y);
y =d * (x^2)*y;

function y=g(x)
y=x^3;

function gr=grad(fnam, x, varargin)
h=0.0000001;

fnam
x
varargin{:}

f1=feval( fnam , x , varargin{:});
f2=feval(fnam,x+h, varargin{:});
gr=(f2-f1)/h
