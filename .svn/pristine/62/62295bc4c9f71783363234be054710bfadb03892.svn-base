function reprocess(settings, infile, outdir)
% reprocess applies the settings to all .nc files under and including
% infile (which may be .nc file or a directory) and writes the output to a
% matching directory structure under outfile.
%
% Inputs:
%   settings    Structure containing settings for process_BASOOP
%   infile      Name of file or directory to reprocess
%   outdir      Directory to write output to [directory of infile].
%
% Author:   Gordon.Keith@csiro.au 20150424

% resolve parameters
% Resolve or ask for parameters
if nargin < 1 
    settings = '';
end

if ischar(settings)
    if isempty(settings)
        [file,pth] = uigetfile('*.mat');
        settings = fullfile(pth,file);
    end
    settings = load(settings);
    if isfield(settings,'settings')
        settings = settings.settings;
    end
end

if nargin < 2 || isempty(infile)
    [file,pth] = uigetfile('*.nc');
    infile = fullfile(pth,file);
end

if nargin < 3 || isempty(outdir)
    outdir = fileparts(infile);
end

% infile is a cell array of netcdf files, redo each one
if iscell(infile)
    for i = 1:length(infile)
        reprocess(settings,infile{i},outdir);
    end
    
% if a directory process all contents of the directory
elseif exist(infile, 'dir') == 7
    if exist(outdir,'dir') == 0
        mkdir(outdir);
    end
    files = dir(infile);
    for i = 1:length(files)
        file = files(i).name;
        if file(1) ~= '.'
            reprocess(settings, fullfile(infile, file), fullfile(outdir, file))
        end
    end
    
% process a netcdf file 
else
    if strcmpi(infile(end-2:end), '.nc')
        settings.read_netcdf = true;
        settings.netcdf_file = infile;
        settings.netcdf = true;
        settings.netcdf_directory = outdir;
        try
            process_BASOOP(settings);
        catch exception
            fprintf('Problem reprocessing %s:\n %s\n', infile, exception.message)
        end
    else
        fprintf('Skipping %s\n', infile)
    end
end

