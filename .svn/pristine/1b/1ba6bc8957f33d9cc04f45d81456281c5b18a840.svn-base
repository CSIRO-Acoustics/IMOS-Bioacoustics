function re_gui(settings, infile, outdir)
% re_gui opens the processing GUI with the settings and netcdf file given
%
% This function will read the list of ev files from the netcdf file and
% create a filelist of them.
%
% Inputs:
%   settings    Structure containing settings for process_BASOOP, 
%               name of a .mat file containing settings or 
%               an empty string will ask the user to select the .mat file
%   infile      Name of netcdf file or directory to reprocess
%   outdir      Location to put output netcdf file
%
% Author:   Gordon.Keith@csiro.au 20150601

% Ensure IMOS toolbox is in path - needed to read netcdf file
imos_path = fileparts(which('imosToolbox'));
if isempty(imos_path)
    imos_path = fullfile(fileparts(mfilename('fullpath')),'IMOS-toolbox');
end

if isempty(strfind(path,imos_path))
    addpath(imos_path);
    addpath(fullfile(imos_path, 'NetCDF'));
    addpath(fullfile(imos_path, 'Parser'));
    addpath(fullfile(imos_path, 'Util'));
    addpath(fullfile(imos_path, 'IMOS'));
    addpath(fullfile(imos_path, 'GUI'));
end
    
% Resolve or ask for parameters
if nargin < 1 
    settings = struct();
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

settings = basoop(settings);

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
        re_gui(settings,infile{i},outdir);
    end
    
% infile is a directory, redo each subdirectory and netcdf file
elseif exist(infile, 'dir') == 7
    if exist(outdir,'dir') == 0
        mkdir(outdir);
    end
    files = dir(infile);
    for i = 1:length(files)
        file = files(i).name;
        if file(1) ~= '.'
            re_gui(settings, fullfile(infile, file), fullfile(outdir, file))
        end
    end
    
% infile is a netcdf file
else
    if strcmpi(infile(end-2:end), '.nc')
        % read netcdf file
        cd(imos_path)
        ncdata = netcdfParse({infile});
        
        % locate ev filenames
        evv = [];
        for i = 1:length(ncdata.variables)
            if strcmp(ncdata.variables{i}.name, 'EV_FILENAME')
                evv = i;
            end
        end
        
        if isempty(evv)
            for i = 1:length(ncdata.dimensions)
                if strcmp(ncdata.dimensions{i}.name, 'EV_FILENAME')
                    evv = i;
                end
                if strcmp(ncdata.dimensions{i}.name, 'CHANNEL')
                    settings.channel = ncdata.dimensions{i}.data;
                end
            end
            if isempty(evv)
                error('Could not find EV_FILENAME in %s', infile);
            end
            evfilelist = ncdata.dimensions{evv}.data;
        else
            evfilelist = ncdata.variables{evv}.data;
        end
        
        % write ev file names to filelist
        [pth,name] = fileparts(infile);
        filelist = fullfile(pth,[name '.ev.txt']);
        fid = fopen(filelist, 'w');
        for i = 1:length(evfilelist)
            evfile = evfilelist{i};
            evfile(evfile == '\') = filesep;
            fprintf(fid,'%s\n', evfile);
        end
        fclose(fid);
    
        % set settings
        evdir = fileparts(evfilelist{1});
        
        settings.alt_ev_files = filelist;
        if ~isempty(evdir)
            settings.echointegration_path = fullfile(evdir, settings.echointegration_directory);
        end
        settings.read_netcdf = true;
        settings.netcdf_file = infile;
        settings.copy_netcdf_metadata = true;
        settings.netcdf = true;
        settings.netcdf_directory = fileparts(outdir);
        
        % show GUI to user 
        try
            GUI(settings);
        catch exception
            fprintf('Problem reprocessing %s:\n %s\n', infile, exception.message)
            keyboard
        end
    else
        fprintf('Skipping %s\n', infile)
    end
end

