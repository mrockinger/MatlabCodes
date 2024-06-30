function u=UnifUnique(N,a,b);
% draws N unique indices in the set [a,b]
% output is a column vector of the indices

% test plausibility of parameters
rge = b-a+1; 
if N>rge
    error('Attempt to generate more unique numbers than possible')
elseif N <= rge
    % if N==rge this corresponds to a random shuffle of the possible numbers
    % the same trick can be nonetheless used for a subsample
    x        = rand(rge,1);
    [sx,idx] = sort(x,1); % really only need the index  
    v=a:1:b; v=v';
    idx=idx(1:N); %take just those indices that are required
    u=v(idx);
end
    
    
end