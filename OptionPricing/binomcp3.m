clc;

% BinomCP3.prg
%
% Estimates binomial tree for a set of values. This allows one to do a comparative static exercise.

CallPut=1; % Positive is call and negative is put

 %S0 K T rho sigma N 
ParmM=[
100 100 60 0.04 0.18 10;
100 100 60 0.04 0.18 20;
110 100 60 0.04 0.18 20;
100 105 60 0.04 0.18 20;
100 100 100 0.04 0.18 20;
100 100 60 0.06 0.18 20;
100 100 60 0.04 0.20 20
];


NNbr=size(ParmM,1);	
CPv=zeros(NNbr,1);  	% prices computed with binomial tree model
CPBSv=zeros(NNbr,1); 	% prices computed with Black-Scholes model

for j=1:NNbr
     
	S0=ParmM(j,1);
	K=ParmM(j,2);
	T=ParmM(j,3)/365;
	rho=ParmM(j,4); 
	rho=log(1+rho);
	sig=ParmM(j,5);
	N=ParmM(j,6);

    fprintf(1, 'S0=%8.4f ,K=%8.4f ,rho=%8.4f ,sig=%8.4f ,T=%8.4f ,N=%8.4f ,CallPut=%2.0f\n',[S0 K rho sig T N CallPut]);

    [StM,CtM,CPrice]=BinomTree(S0,K,rho,sig,T,N,CallPut);
    C=BScall(S0,K,sig,T,rho);

	CPv(j)=CPrice;
	CPBSv(j)=C;

end;

format short;
fid=fopen('d:/finbox/option/out.out','w');
%fid=1;
ResM=[ParmM CPv CPBSv];
fmt0='%12.4f';
fmt=Kron(ones(1,size(ResM,2)),fmt0); 
fmt=[fmt '\n'];
fprintf(fid,fmt,ResM);
if fid>1;
    fid=fclose(fid);
end   