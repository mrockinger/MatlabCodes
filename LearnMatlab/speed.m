function speed()
% illustrates how to write a program and how not to do it

% Do !

T=10000;
tic;
Res=zeros(T,1);
for j=1:2
for i=1:T
    Res(i,1)=i+j;
    if i>2
        break
    end
end
end
100000*toc

% Do NOT!

tic
Res=[];
for i=1:T
    Res=[Res; randn(1,1)];
end
100000*toc

tic
Res=randn(T,1);
100000*toc