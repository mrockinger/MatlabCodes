function TraceSimplex()
% traces a simplex to illustrate how a grid may get generalized into two
% dimensions
clc

a1=0:0.1:1; a1=a1'
Na1=rows(a1)
a1r=[];
a2r=[];
a3r=[];

%hold on;
for i=1:Na1
    a2=0:0.1:1-a1(i); a2=a2'
    a1r=[a1r; a1(i).*ones(size(a2))];
    a2r=[a2r; a2];
end    
a3r=1-a1r-a2r;
[a2r a1r  a3r]
plot3(a1r,a2r,a3r,'ks')    
hold on;
t=0:0.1:1; t=t'; Nt=rows(t);
 plot3(zeros(Nt,1),t,1-t,'k');
 plot3(t,zeros(Nt,1),1-t,'k');
 plot3(t,1-t,zeros(Nt,1),'k');
hold off;
title('Three dimensional grid descirbing the weighting of the mixture densities')
xlabel('\alpha_1');
ylabel('\alpha_2');
zlabel('\alpha_3');
