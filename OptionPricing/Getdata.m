% loads data and checks that it is OK

X=xlsread('StkIdx.xls'); 
[r,c]=size(X);
fprintf('there are %4.0f rows and %4.0f columns\n', [r c]);

% now verify quality of data
fmt=['%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f \n'];
fprintf(fmt,X(1,:));
fprintf(fmt,X(2,:));
fprintf(fmt,X(3,:));
fprintf('\n');
fprintf(fmt,X(r-2,:));
fprintf(fmt,X(r-1,:));
fprintf(fmt,X(r,:));

