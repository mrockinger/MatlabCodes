function [V,CV]=get_Var_emp(x,Conf_lev,LR);
% computes the VaR, V, and the CVaR, CV, for a given matrix
% il LR<=0 look at left tail. Else takes right tail.
[T,K]=size(x);

idx=floor(Conf_lev*T);
if idx==0
    error('Error in get_Var_emp: Conf_level too small. Index=0 ')
end

V  = zeros(K,1);
CV = zeros(K,1);

if LR>=0
    x=-x;
end

sx = sort(x);
V  = sx(idx,:);
CV = mean(sx(1:idx,:));

if LR>=0
    V=-V; CV=-CV;
end