function reprocessev(settings, infile, outdir)
% reprocessev applies the settings to all .nc files under and including
% infile (which may be .nc file or a directory) and writes the output to a
% matching directory structure under outfile.
%
% It also reads the list of evfiles from the original netcdf file and uses
% that list for any processing involving evfiles.
% 
% If the create ev files option is selected then the template specified will be
% used to create new ev files with the new template and the original raw
% files with all regions and lines copied from the original ev files.
%
% If a new netcdf file is written then the metadata from the original
% netcdf file will be copied.
%
% Inputs:
%   settings    Structure containing settings for process_BASOOP
%   infile      Name of file or directory to reprocess. 
%   outdir      Directory to write ouput netcdf files. [settings.netcdf_directory]
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
% resolve parameters
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
    if isfield(settings, 'netcdf_directory')         
        outdir = settings.netcdf_directory;        
    else
        outdir = fileparts(infile);                
    end
end


% infile is a cell array of netcdf files, redo each one
if iscell(infile)
    for i = 1:length(infile)
        reprocessev(settings,infile{i},outdir);
    end  
% if a directory process all contents of the directory
elseif exist(infile, 'dir') == 7
    if exist(outdir,'dir') == 0
        mkdir(outdir);
    end
    
    files = dir(infile); 
    for i = 1:length(files)
        file = files(i).name;
        if file(1) ~= '.' && isequal(files(i).name(end-1:end),'nc'); % changed to filter for nc files   
            
            % reprocessev(settings, fullfile(infile, file),
            % fullfile(outdir, file)) % blatted this out as we want the
            % output directory not the directory + nc file (TR 07/08/2017)
            reprocessev(settings, fullfile(infile, file), fullfile(outdir))
        end
    end
    
% process a netcdf file 
else    
    if strcmpi(infile(end-2:end), '.nc')
        % read netcdf file
        fprintf('%s: Reading file %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), infile);
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
                    % settings.channel = ncdata.dimensions{i}.data;
                    % TER changed to convert channels from vertical vector
                    % to horizontal
                    settings.channel = ncdata.dimensions{i}.data';
                end
            end
            if isempty(evv)
                warning('Could not find EV_FILENAME in %s', infile);
                evfilelist = {};
            elseif ~settings.use_alt_ev_files    & settings.read_netcdf            
                evfilelist = ncdata.dimensions{evv}.data;
                if ~ exist(evfilelist{1})
                    warning('Could not find EV_FILENAME in %s\n',evfilelist{1})
                    warning('Most likely path embedded in netcdf has changed');                
                end
                % ok, need to provide a new path to the ev file as has
                % probably changed since processing. Tim  Ryan 20/08/2017
                % netcdf file has evfile list, but the evfile list is pointing to
                % defunct locations. there is a good chance the evfiles will be
                % located alongside the ncdfile 
                for i=1:size(evfilelist,1)
                    if ~exist(evfilelist{i})% ev file not found. remap the path
                        sep = strfind(infile,'\');
                        nc_root_dir = infile(1:sep(end));
                        evworksheet_dir = [nc_root_dir 'Echoview_worksheets'];
                        [evdir,evfilename,evext] = fileparts(evfilelist{i});
                        evfile = [evworksheet_dir '\' evfilename evext];
                        evfilelist{i} = evfile;                            
                        if ~ exist(evfile)
                            fprintf('Ev file not deduced from netCDF information\n')
                            fprintf('Provide this by pointing to a text file in the GUI where it says\n')
                            fprintf('"Use ev files from" ... \n')
                            fprintf('This is a text file that provides the path for each ev file that you wish to process\n')
                        end
                    end
                end                
% % % %             elseif isempty(ls(settings.alt_ev_files)) & ~ exist(evfile)
% % % %                     warning(sprintf(['\n=============================================================\n' ...
% % % %                                     'Alternative ev file text file does not exist.\n' ...
% % % %                                     'Will try to figure out location of ev files from netcdf info\n',...
% % % %                                     'otherwise go set this in settings.alt_ev_files (use the GUI or do this as text entry)\n'...
% % % %                                     '=================================================================\n']))                                              
            elseif settings.use_alt_ev_files    & ~settings.read_netcdf                       
                idx =1;                
                fid1 = fopen(settings.alt_ev_files); 
                if isequal(fid1, -1)
                    fprintf('alt_ev_file not found. program will crash. go specify an alt_ev_file\n')
                else
                    while ~ feof(fid1)
                        ln = fgetl(fid1);                    
                        evfilelist{idx} = ln;

                        if isempty(ls(ln))
                            fprintf('Ev file not deduced from netCDF information\n')
                            fprintf('Nor is it found in the text file provided in the \n')
                            fprintf('"Use ev files from" ... dialog\n')
                            fprintf('This is a text file that provides the path for each ev file that you wish to process\n')
                            fprintf('Go make sure you have the right text file that points to ev files that you wish to reprocess\n');
                        end
                        idx = idx+1;
                    end     
                    fclose(fid1);
                end
            end
        
        else                     
            evfilelist = ncdata.variables{evv}.data;
        end             
        % Some older data sets had the filenames transposed
        % e.g. 'QQQ' not 'Q:\'
        % detect and fix this problem
       
        if ~isempty(evfilelist) && evfilelist{1}(1) == evfilelist{1}(3)
            flist = cell2mat(evfilelist);
            filelist=reshape(flist', size(flist));
            for i = 1:length(evfilelist)
                evfilelist{i} = filelist(i,:); %#ok<AGROW>
            end
        end        
        % write ev file names to filelist
        [pth,name] = fileparts(infile);
        filelist = fullfile(pth,[name '.ev.txt']);
        fid = fopen(filelist, 'w');
        for i = 1:length(evfilelist)
            evfile = evfilelist{i};
 %           evfile(evfile == '\') = filesep;
            fprintf(fid,'%s\n', evfile);
        end
        fclose(fid);
            
        % set settings                
        % ----  modifications to this section by Tim Ryan 09/08/2017 ----
        %
        if isempty(evfilelist) || isempty(ls(evfilelist{1}))  % ev files list embedded in the netcdf cannot be found. possibly path has changed since processing. look to use the file list specified in the settings.alt_ev_files 
            if ~isempty(ls(settings.alt_ev_files)) 
            % do nothing - settings.alt_ev_files picked up from the
            % settings.alt_ev_files in the settings structure (possibly as
            % per the GUI entry)
            filelist = settings.alt_ev_files
            else
               fprintf('hmmm, ev files list embedded in the netcdf file cannot be found, nor can evfiles listing in the file in the settings.alt_ev_files value be found\n');
               fprintf('solution: nothing you can do about the list of files in the netcdf, but you can provide a list to settings.alt_ev_files that does have the correct path for the ev files\n');
               fprintf('program will now crash - \n');               
            end
        else % ev files embedded in the netcdf are found - program can continue :) 
            settings.alt_ev_files = filelist;
        end
        % ---------------- end modifications -----------------------
        settings.overwrite_ev_files = true;
        settings.import_ev = [];
        settings.netcdf_file = infile;
        if settings.netcdf
            settings.read_netcdf = true;
            settings.copy_netcdf_metadata = true;
            settings.update_platform = [];            
            settings.netcdf_directory = outdir;
        end        
        % if create .ev files then copy regions from existing .ev files
        if settings.create_ev_files
            try
                EvApp = actxserver(settings.EvApp);
            catch e
                EvApp = [];
                fprintf('Couldn''t create ActiveX server for echoview\n%s\n', e.message);
            end
            
            datalist = {{}};
            
            for f = 1:length(evfilelist)
                fprintf('%s: Opening file %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), evfilelist{f});
                EvFile = EvApp.OpenFile(evfilelist{f});
                if ~isempty(EvFile)
                    [pth,file,~] = fileparts(evfilelist{f});
                   
                    for s = 1:length(settings.filesets)
                        dataset = EvFile.Filesets.FindByName(settings.filesets{s});
                        if ~isempty(dataset)
                            for j=1:dataset.DataFiles.Count
                                filename = dataset.DataFiles.Item(j-1).FileName;
                                if length(datalist) < s || isempty(datalist{s}) || ~strcmp(datalist{s}{end}, filename)
                                    datalist{s}{end+1} = filename;
                                end
                            end
                        end
                    end
                    
                    % If we don't find the filesets named in
                    % settings.filesets then use the first fileset
                    if isempty(datalist{1})
                       dataset =  EvFile.Filesets.Item(0); 
                       settings.filesets{1} = dataset.Name;
                       for j=1:dataset.DataFiles.Count
                           filename = dataset.DataFiles.Item(j-1).FileName;
                           if isempty(datalist{1}) || ~strcmp(datalist{1}{end}, filename)
                               datalist{1}{end+1} = filename;
                           end
                       end
                    end
                    
                    dataset = EvFile.Filesets.FindByName('Transit_GPS');
                    if ~isempty(dataset) && dataset.DataFiles.Count > 0
                        settings.transit_gps_file = dataset.DataFiles.Item(0).FileName;
                    end
            
                    dataset = EvFile.Filesets.FindByName('Roll');
                    if ~isempty(dataset) && dataset.DataFiles.Count > 0
                        settings.transit_roll_file = dataset.DataFiles.Item(0).FileName;
                        settings.include_roll = true;
                    end
            
                    dataset  = EvFile.Filesets.FindByName('Pitch');
                    if ~isempty(dataset) && dataset.DataFiles.Count > 0
                        settings.transit_pitch_file = dataset.DataFiles.Item(0).FileName;
                        settings.include_pitch = true;
                    end
                    
                    exports = {};
                    fprintf('%s: Exporting regions\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
                    evr = fullfile(pth, [file '.evr']);
                    if EvFile.Regions.ExportDefinitionsAll(evr);
                        exports{end+1} = evr;           %#ok<AGROW>
                    else
                        warning('Failed to export regions from %s\n', evfilelist{f});
                    end
                    fprintf('%s: Exporting lines\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
                    for l =0:EvFile.Lines.Count-1
                        linename=EvFile.Lines.Item(l).Name;
                        if isempty(find(linename == ':',1))
                            evl = fullfile(pth, [file '.' linename '.evl']);
                            if EvFile.Lines.Item(l).Export(evl)
                                exports{end+1} = evl;       %#ok<AGROW>
                            end
                        end
                    end
                    
                    EvFile.Close();
                else
                    keyboard
                    error('Unable to open EV file %s', evfilelist{f})
                end
                
                settings.import_ev{f} = exports;
            end
            
            EvApp.Quit();
            
            for f = 1:length(datalist)
                listfile = fullfile(fileparts(infile),sprintf('datalist%d.txt', f));
                fid = fopen(listfile, 'w');
                for j = 1:length(datalist{f})
                    fprintf(fid,'%s\n', datalist{f}{j});
                end
                fclose(fid);
               
                try
                    settings.transit_data_files{f} = listfile;
                    % crashing - don't know why. Tim Ryan 24/08/2017
                catch
                    settings.transit_data_files = listfile; 
                end
            end
            
        end        
        % process
        fprintf('%s: Starting reprocess for %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), infile);
        try                     
            process_BASOOP(settings);
        catch exception
            fprintf('%s: Problem reprocessing %s:\n %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
                infile, exception.message)
            display(exception.stack(1))
            keyboard
        end
    else
        fprintf('Skipping %s\n', infile)
    end
end

