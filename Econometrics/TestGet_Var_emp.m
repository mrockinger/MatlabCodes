% Test_Get_Var_Emp

x=randn(10,3);
x=sort(x) % to better see later on what goes on
Conf_lev=0.05;
Conf_lev=0.25;
LR=-1; % left tail

[V,CV]=get_Var_emp(x,Conf_lev,LR);

disp('VaR level');
V
disp('CVaR');
CV

[V,CV]=get_Var_emp(x,Conf_lev,-LR);

disp('VaR level');
V
disp('CVaR');
CV

[V,CV]=get_Var_emp(x,0.0001,-LR);
