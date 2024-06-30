    function learning3()
    % learning3.m
    % a first matlab function
    % loads some ascii data, 
    % creates with that data a Matlab dataset. Loads data and displays it
    %
    clc;
    fid = fopen('d:\\finbox\\my2.txt','r')
    [A,T] = fscanf(fid,'%f %f %f'); % no missing values.
    fclose(fid);
    
    A
    
    fid = fopen('d:\\finbox\\my2.txt','r')
    [A,T] = fscanf(fid,'%f %f %f',[3 6]); % no missing values 
    fclose(fid);
    
    A
    % !!! A is transposed
    A=A'
    
    format short;
    A(1:3)
    A(end-2:end)
    disp('T');
    T
    
    x=reshape(A,3,T/3)' %create 3 rows then transpose Each unit has T/3 elements
    x
    
    % now writes data into a new file
    disp('Save data into file');
    save 'd:\finbox\my3' x;
    whos

    clear;
    disp('whos after clear')
    whos
    
    tic
    load 'd:\finbox\my3' x
    x
toc
return

    % save as an ascii file    
    fid = fopen('d:\\finbox\\my3.txt','w');
    fprintf(fid,'%5.0f %8.3f %8.3f \n',x); %don't forget the \n
    % notice \n is readable under word but not under wordpad...
    fclose(fid);
   
    
    