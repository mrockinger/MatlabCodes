function learning1()
% learning Matlab
echo on
A=[ 1 2 3; 4 5 6] % creates a Matrix

b=[88; 99]; % notice the difference due to ;

disp('a column vector');
b

fprintf('another way to print b but using formatted output. b=%12.4f\n',b);


disp('what is the precision of Matlab');
fprintf('pi=%25.20f',pi);

A(:,1) % the first column of A
A(2,:) % the second row of A
A(1,3) % an element out of A
A(2,2:3)=[-2 -3]; % replace blocks in A
A
% A(2,1 3)=[5 7]
% return

B=[1 7 3 ; 5 6 7] % another matrix
B' %transpose of B

C=B'*A %yet another matrix
5*C % multiply a matrix by a scalar


D=[A;B] % vertical concatenation (combine vertically A and B)
D=[A B] % horizontal concatenation

A.*B % Haddamar product (element by element product)

E=[] % an empty matrix
E=[E;A;A] %concatenate matrices (useful to store results in E)

A+B %sum of matrices
A-B % difference of matrices

a=A(1,:)
b=B(:,2)

[r,c]=size(A) %rows and columns of A
size(A,1) % row dimension of A
size(A,2) % column dimension of A


% notice: a.*b would return a mistake
kron(a,b)

clear a
%just typing a  returns an error

A=1:3 % creates a sequence for A. Notice A is a row!!!

A(end)
A(end-1)


B='I am a string of characters'
C='; lots of them'
D=[B C] % horizontal concatenation

S=char(128) % the Euro
D=[D ' ' S] 

char(65)

E=strvcat(B,'Yes') % vertical concatenation.  Notice padding
disp(E)

E=strcat(B,' No') % equivalent to [B ' No']. Emphasis on string operation

double('A')  % to go from char to integer
double('a')


x=randn(3,1) % generates normal random numbers
str1 = num2str(min(x)) % a first string
str2 = num2str(max(x)) % another one
disp(str1); %only strings can get displayed with disp
['from smallest ' str1 ' to largest ' str2] 

s=num2str(min(x),12) % more decimals
['now display 12 decimals ',s]

str = '12.3389e-1';
val = str2num(str) % goes from string to real
val+12

% do not confuse with cells. These are structures.
A={[1; 2] 'I scream for ice cream' }
A{1,1}
A{1,2}
return
