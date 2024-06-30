function WhileEx()
% WhileEx illustrates the use of the while loop


i=1;
while i<10
    [i i^2]
    i=i+1;
end



Res=[];
A=zeros(10,1);
for x=1:0.1:2
     Res=[Res; [x x^2]]; % only ok for short series else pre-declare size
     A(floor(x*10),:)=x
end
Res
disp('A');
A