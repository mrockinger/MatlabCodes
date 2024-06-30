function hrulefill(n, fid)

if nargin < 1
    n = 80;
end
if nargin < 2
    fid = 1;
end
fprintf(fid, '%s\n', repmat('-', 1, n));