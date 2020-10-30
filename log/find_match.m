function matches = find_match(match, directory)
% FIND locates all the .mat files in the directory which have the fields
% and values found in match.
%
% written to scan through the IMOS_BASOOP log directory to find usable test
% cases

if nargin < 2 || isempty(directory)
    directory = '.';
end

matches={};
fields = fieldnames(match);

files = dir(fullfile(directory, '*.mat'));

for i = 1 : length(files)
    try
    test = load(fullfile(directory,files(i).name));
    catch e
        e.message
    end
    ok = true;
    for field = fields';        
        if ~isfield(test, field{1}) 
            ok = false;
            break;
        end
        if test.(field{1}) == match.(field{1})
        else
            ok = false;
            break;
        end
    end
    if ~ok
        continue;
    end
    matches{end+1} = files(i).name;
end

    