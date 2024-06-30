    function fcall()
    % illustrates how to nest functions and local/globalness of variables
    global c
    a=2
    b=3
    c=10
    
    [x,y]=do_Comput(a,b) %here involves two elements that are returned
   
    
    do_more_Comput(a,b)
    
    function [f,g]=do_Comput(a,b);
    global c
    f=a^2+c; 
    g=a*b*c;
    
    