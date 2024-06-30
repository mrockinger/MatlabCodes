% the usual mathematical functions
clc
echo on
A=[1 2 3; 1 2 3]
rank(A) %rank of A
B=[1 2; 3 4]
rank(B)
det(B) %determinant
trace(B) %trace
[V,D] = eig(B) %eigen vectors and eigen values
[B*V(:,1) D(1,1)*V(:,1)] %verification that B v_1= lda_1 v_1

inv(B)
c=[3; 4]
inv(B)*c % one way to compute
% a much better way to go (for numerical precision)
B\c
E=eye(2)
% if one needs to invert B use (faster + more precise)
B\E

% condition number (ratio of largest to smallest eigenvalue)
% if larger than 30'000 be careful when invering matrix
cond(B)
C=[1 2; 0.9999 2.0001]
cond(C)

%Choleski decomposition matrix must be positive definite
B=[1 0.1; 0.1 3]
X=chol(B)
X'*X %verification
