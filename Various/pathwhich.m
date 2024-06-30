% recurse through path and check for overlapping functions
%initscript

% note: the precedence order is determined from the point of view of the current working directory

% some obvious files can be excluded:
excludefiles = {'contents.m', 'make_html.m'};

p = path;
% If I am only interested in my self-defined path, that is only what's preceding mathworks
p = p(1 : strfind(path, matlabroot) - 1);

[thisdir, restdirs] = strtok(p, pathsep);
overlaps = cell(0);
while ~isempty(restdirs)
    whatfiles  = what(thisdir);
    for m = 1 :length(whatfiles.m)
        if isempty(strmatch(lower(whatfiles.m{m}), lower(excludefiles)))
            c = which(whatfiles.m{m}, '-all');
            % todo: prune @ overloads
            if length(c) > 1 && length(c) ~= length(strmatch(matlabroot, c)) % final clause: skip overloads caused by mathworks
                if length(overlaps) == 0 || isempty(strmatch(whatfiles.m{m}, overlaps, 'exact')) ...
                        overlaps{end+1} = whatfiles.m{m}; % no double counting
                    hrulefill;
                    fprintf('Overloaded function %s. Instances are:\n', whatfiles.m{m});
                    fprintf('%s (active)\n', c{1});
                    fprintf('%s \n', c{2:end});
                    hrulefill;
                end
            end
        end
    end
    [thisdir, restdirs] = strtok(restdirs, pathsep);
end

% Final summary
%hrulefill
fprintf('Found %d overshadowed functions.\n', length(overlaps))
fprintf('There names are:\n')
[junk, ndx]     = sort(lower(overlaps));
fprintf('%s\n', overlaps{ndx})

