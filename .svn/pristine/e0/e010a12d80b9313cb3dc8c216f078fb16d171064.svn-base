function ev_file_list = create_ev_files(file_sets, EvApp, control, progress)
% Creates EchoView worksheets for the the specified data files.
%
% Inputs:
%   file_sets   cell array of cell arrays of ES60 (or EK500) file names, an
%               .ev file is created for each row of the outer array 
%               containing all the files in the inner array.
%               The number of columns in file_sets equals the length of
%               control.filesets.
%   EvApp       a handle to a COM object of an EchoviewApplication.
%   control     structure containing the following fields:
%       echoview_file_location  directory to put .ev files
%       template                EchoView template used to create worksheets
%       transit_gps_file        name of file containing GPS for the transit
%       detect_seafloor         flag for seafloor detection
%       detect_upper_DSL        flag for deep scattering layer detection
%       detect_fixed            flag for fixed depth layer
%       fixed_layer             depth of fixed depth layer
%       upper_DSL_line          name of upper DSL or fixed line
%       frequency               transducer frequency
%       base_variable_name      'Sv_%gkHz'
%       resample_variable_name  'Sv %g kHz resample for DSL convolution'
%       filesets                cell array of Fileset names in the template
%                               to be populated with files from file_sets
%                               for example {'Vessel_sv_data'}
%       time_offset             difference between raw file times and UTC
%       include_pitch           flag to include pitch.csv file
%       transit_pitch_file      name of file containing pitch for the
%                               transit, only required if include_pitch
%       include_roll            flag to include roll.csv file
%       transit_roll_file       name of file containing roll for the
%                               transit, only required if include_roll
%       read_ecs                flag to use ecs file
%       calibration_file        .ecs file to use
%   progress    Optional handle to a progress function progress(i,n,file)
%
% Ouputs:
%   ev_file_list name of file containing list of worksheets created.
%
% Author:
%       Tim Ryan <tim.ryan@csiro.au>
%       Gordon Keith <gordon.keith@csiro.au>
%       2011-10-18

    if nargin < 4
        progress = [];
    end
    
    % set of lines to create in the ev file.
    fixed_lines = struct(...
        'name', { 'exclude firepulse', 'Upper Line', 'Lower Line', '120kHz Lower Line', 'Upper DSL line'}, ...
        'depth', { 10, 0, 250, 80, 500});
    
    % if echoview connection is not provided run echoview
    if isempty(EvApp)
        EvApp = actxserver('EchoviewCom.EvApplication');
    end
    last_logged = EvApp.GetLastLogMessage();
    EvApp.Minimize; % Haris 25 July 2019 - minimize Echoview window opened
    
    ev_file_list = ['Ev_files_' simrad_date_string(file_sets{1,1}{1}) '.txt'];
    ev_file_list = fullfile(control.echoview_file_location, ev_file_list);
    
    % copy the template file to the processed data directory just once
    % and name with the date that the data was processed
    template = fullfile(control.echoview_file_location, ...
        ['Template_D' datestr(now, 'yyyymmdd_THHMMSS') '.ev' ]);
    copyfile(control.template, template);
    
    fid = fopen(ev_file_list, 'w');
    
    for i=1:size(file_sets,1);
        if ~isempty(progress)
%             progress(i,size(file_sets,1),file_sets{i,1}{1}); %#ok<NOEFF>
            progress(file_sets{i,1}{1}); % now displaying rawfile name and message is clear -Haris 10 July 2019
        end
        
        % ---- create the ev file name based on the name of the first file
        % in each of the file_sets
        evfilename = [simrad_date_string(file_sets{i,1}{1}) '.ev' ];        
        evfilename = fullfile(control.echoview_file_location, evfilename);
        
        if ~control.overwrite_ev_files && exist(evfilename, 'file') == 2
            fprintf(fid, '%s\n', evfilename); 
            continue
        end
        
        % ---- create the EV file based on named template file ----
        EvFile = EvApp.NewFile(template);
        if EvFile.SaveAs(evfilename);
            fprintf('                     created Ev file %s\n', evfilename);
        else
            last_logged = print_new_log(EvApp, last_logged);
            fprintf('Failed to create ev file %s \n', evfilename);
        end
        fprintf(fid, '%s\n', evfilename); % collate list of created files in a text file for later reference. 
        pause(1);
        
        % first add the transit gps file so that the template can reference this
        % before adding the raw files     
        gps = EvFile.Filesets.FindByName('Transit_GPS');
        if control.time_block == 0 && (isempty(gps) || exist(control.transit_gps_file, 'file') ~=2)
            % gps.csv is optional if time_block is 0
        else
            if isempty(gps)
                print_new_log(EvApp, last_logged);
                error('Transit_GPS fileset not found in template');
            end
            if ~gps.DataFiles.Add(control.transit_gps_file);
                last_logged = print_new_log(EvApp, last_logged);
                fprintf('GPS file %s cannot be found. Check that this file has been created and is in the correct location\n', control.transit_gps_file)
            end
        end
        gps.TimeOffset = control.time_offset(1);
        
        roll = EvFile.Filesets.FindByName('Roll');
        if control.include_roll
            if isempty(roll)
                print_new_log(EvApp, last_logged);
                error('Roll fileset not found in template');
            end
            if ~roll.DataFiles.Add(control.transit_roll_file);
                last_logged = print_new_log(EvApp, last_logged);
                fprintf('Roll file %s cannot be found. Check that this file has been created and is in the correct location\n', control.transit_roll_file)
            end
            roll.TimeOffset = control.time_offset(1);
        end
        
        pitch = EvFile.Filesets.FindByName('Pitch');
        if control.include_pitch
            if isempty(pitch)
               print_new_log(EvApp, last_logged);
               error('Pitch fileset not found in template');
            end           
            if ~pitch.DataFiles.Add(control.transit_pitch_file);
                last_logged = print_new_log(EvApp, last_logged);
                fprintf('Pitch file %s cannot be found. Check that this file has been created and is in the correct location\n', control.transit_pitch_file)
            end
            pitch.TimeOffset = control.time_offset(1);
        end
        
        % Use dummy pitch and roll if none provided
        if ~isempty(roll)
            if roll.DataFiles.Count == 0
               rollfile = 'empty.roll.csv';
               if isfield(control, 'default_roll_file')
                   rollfile = control.default_roll_file;
               end
               if exist(rollfile, 'file') ~= 2
                   rollfile = fullfile(fileparts(mfilename('fullpath')),rollfile);
               end
               roll.DataFiles.Add(rollfile);
            end    
        end
        
        if ~isempty(pitch)
            if pitch.DataFiles.Count == 0
               pitchfile = 'empty.pitch.csv';
               if isfield(control, 'default_pitch_file')
                   pitchfile = control.default_pitch_file;
               end
               if exist(pitchfile, 'file') ~= 2
                   pitchfile = fullfile(fileparts(mfilename('fullpath')),pitchfile);
               end
               pitch.DataFiles.Add(pitchfile);
            end    
        end
        
        % ---- add listed set of raw files to the fileset
        for f = 1:size(file_sets,2)
            dataset = EvFile.Filesets.FindByName(control.filesets{f});
            if isempty(dataset) && f == 1           % ensure backward compatibility, may not be needed
                dataset =  EvFile.Filesets.Item(0);
            end
            if ~isempty(dataset)
%                 fprintf('%s: adding %d files to %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
%                     length(file_sets{i,f}), dataset.Name);
                fprintf('                     adding %d files to %s\n', length(file_sets{i,f}), dataset.Name);
                if f < length(control.time_offset)
                    dataset.TimeOffset = control.time_offset(f);
                else
                    dataset.TimeOffset = control.time_offset(1);
                end
                if control.read_ecs
                    [~, ecs_filename, ecs_ext] = fileparts(control.calibration_file);
                    ecs_file = fullfile(control.echoview_file_location, [ecs_filename ecs_ext]);
                    copyfile(control.calibration_file, ecs_file);
                    dataset.SetCalibrationFile(ecs_file);
                end
                for j=1:length(file_sets{i,f})
                    if ~dataset.DataFiles.Add(file_sets{i,f}{j})
                        fprintf('Unable to add file %s\nto %s\n', file_sets{i,f}{j}, control.filesets{f});
                    end
                end
            end
        end
        
        last_logged = EvApp.GetLastLogMessage();
        
        % ---- detect seafloor 
        if control.detect_seafloor
            % first set the line pick algorithm and its parameter values
            EvFile.Properties.LinePick.Algorithm = 'eLinePickBestBottomCandidate';
            EvFile.Properties.LinePick.StartDepth = 10;
            EvFile.Properties.LinePick.StopDepth = 1500;
            EvFile.Properties.LinePick.MinSv = -50;
            EvFile.Properties.LinePick.UseBackstep = 1; % i.e true
            EvFile.Properties.LinePick.DiscriminationLevel = -50.0;
            EvFile.Properties.LinePick.BackstepRange = -30; % don't mess around, get the acoustic bottom well away from the seafloor.
            EvFile.Properties.LinePick.PeakThreshold = -35.0;
            EvFile.Properties.LinePick.MaxDropouts = 20;
            EvFile.Properties.LinePick.WindowRadius = 8;
            EvFile.Properties.LinePick.MinPeakAsymmetry = -1.0;   
            fprintf('                     picking acoustic bottom\n');
            base_variable = sprintf(control.base_variable_name, control.channel{1});
            % at some point need some work here to allow frequency specfic
            % acoustic bottom to be defined. Note acoustic bottom to be
            % defined for all channels. So a for loop may be needed for
            % defining acoustic bottom_38kHz, acoustic bottom_120kHz etc.
            % This is what defined in the template                                    
            pick_line(EvFile, base_variable, 'acoustic bottom');    
            last_logged = print_new_log(EvApp, last_logged);
        end
        
        % ----- detect deep scattering layer
        if control.detect_upper_DSL
            fprintf('Picking upper deep scattering layer (DSL) line \n');
            EvFile.Properties.LinePick.Algorithm = 'eLinePickMaxSv';
            EvFile.Properties.LinePick.StartDepth = 250;
            EvFile.Properties.LinePick.StopDepth = 750;
            EvFile.Properties.LinePick.MinSv = -70;
            EvFile.Properties.LinePick.UseBackstep = 1; 
            EvFile.Properties.LinePick.BackstepRange = -50;
            resample_variable = sprintf(control.resample_variable_name, control.channel{1});
            pick_line(EvFile, resample_variable, control.upper_DSL_line)
            last_logged = print_new_log(EvApp, last_logged);
        end          
        
        % ----- fixed line as basis for deep scattering layer
        if control.detect_fixed
            FixedLine = EvFile.Lines.CreateFixedDepth(control.fixed_layer);
            Line =  EvFile.Lines.FindByName(control.upper_DSL_line);
            if isempty(Line)
                warning('CREATEEV:UDSL', 'Line "%s" not found in template', control.upper_DSL_line );
                FixedLine.Name = control.upper_DSL_line;
                EvFile.Lines.Delete(FixedLine);
            else
                Line = Line.AsLineEditable;
                if isempty(Line)
                    error('Line %s is not editable',control.upper_DSL_line);
                end
                Line.OverwriteWith(FixedLine); 
                EvFile.Lines.Delete(FixedLine);
            end
        end
        
        % ------ Create standard fixed lines 
        
        % HK 26/11/2019- adding 'if control.fixed_lines' (default to true
        % in basoop.m).Lines are frequency specific in new fast processing
        % template. Work is needed to make it frequency specific. This
        % setting will be turned off in fast processing to avoid series of
        % warning messages. The template correctly captures all the lines.
        
        if control.fixed_lines 
            for f = 1: length(fixed_lines)
                FixedLine = EvFile.Lines.CreateFixedDepth(fixed_lines(f).depth);
                Line =  EvFile.Lines.FindByName(fixed_lines(f).name);
                if isempty(Line)
                    warning('CREATEEV:FIXED', 'Line "%s" not found in template', fixed_lines(f).name);
                    FixedLine.Name = fixed_lines(f).name;
                    EvFile.Lines.Delete(FixedLine);
                else
                    Line = Line.AsLineEditable;
                    if isempty(Line)
                        error('Line %s is not editable',fixed_lines(f).name);
                    end
                    Line.OverwriteWith(FixedLine); 
                    EvFile.Lines.Delete(FixedLine);
                end
            end
        end
        
        if isfield(control, 'import_ev')
            if iscell(control.import_ev)
                if i <= length(control.import_ev)
                    import_files = control.import_ev{i};
                else
                    warning('CREATEEV:MISSING_IMPORT', 'Creating %d EV files but only %d to import.\nCheck time block is correct', size(file_sets,1), length(control.import_ev));
                end
            else
                import_files = {control.import_ev};
            end
            if ~iscell(import_files) 
                import_files = {import_files};
            end
            for r = 1:length(import_files)
                [~,fname,ext] = fileparts(import_files{r});
                fprintf('%s: Import: %s%s\n',  datestr(now, 'yyyy-mm-dd HH:MM:SS'), fname, ext);
                if strcmp(ext, '.evl')
                    [~,~,lname] = fileparts(fname);
                    if isempty(lname)
                        error('Unable to find line name in "%s"', import_files{r})
                    else
                        lname(1) = []; % drop '.'
                    end
                    line =  EvFile.Lines.FindByName(lname);
                    nlines = EvFile.Lines.Count;
                    if EvFile.Import(import_files{r});
                        if nlines < EvFile.Lines.Count
                            newLine = EvFile.Lines.Item(nlines);
                            if ~isempty(line)
                                line = line.AsLineEditable();
                                if ~isempty(line)
                                    if line.OverwriteWith(newLine);
                                        fprintf('%s: line "%s" successfully updated\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), lname);
                                        if ~EvFile.Lines.Delete(newLine)
                                            warning('CREATEEV:DELETETLINE', 'Deleting of line failed');
                                        end
                                    else
                                        warning('CREATEEV:OVERWRITEFAIL', 'Overwriting line "%s" failed', lname);
                                    end
                                else
                                    warning('CREATEEV:LINENOTEDITABLE', 'Line "%s" not editable (virtual)', lname);
                                end
                            else
                                line.Name = lname;
                            end
                        else
                            warning('CREATEEV:IMPORT', 'Cannot find new line "%s"', lname);
                        end
                    else
                        warning('CREATEEV:IMPORT', 'Import failed for line "%s" \n  file %s', lname, import_files{r});
                    end
                                       
                else
                    if ~EvFile.Import(import_files{r});
                        warning('CREATEEV:IMPORT', 'Import failed for file %s', import_files{r});
                    end
                end
            end
        end
       
        
        EvFile.Save;        
        try
            EvFile.Close;          
        catch exception
            % workaround for echoview 7 intermittant crash on close bug
            warning('CREATEEV:EV7_CLOSE', 'Echoview close problem %s', exception.message)
            EvApp = actxserver('EchoviewCom.EvApplication');
        end
%         fprintf('%.1f percentage complete\n', 100*i/size(file_sets,1));
        fprintf('%s: %.1f percentage complete\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'),100*i/size(file_sets,1));
    end
    fclose(fid); 
end

function pick_line(EvFile,VarName, Line_Name)
    % assumes ev_file has been created using the appropriate template 
    Var = EvFile.Variables.FindByName(VarName); 
    if isempty(Var)
        fprintf('Variable not found\n'); 
    else
        VarAc = Var.AsVariableAcoustic;    
        Lines = EvFile.Lines();
        Line = Lines.CreateLinePick(VarAc, 1);
        
        OldLine = Lines.FindByName(Line_Name);
        if isempty(OldLine)
            warning('CREATEEV:PICK','Line "%s" not found in template', Line_Name);
            Line.Name = Line_Name;
            Lines.Delete(Line);
        else
            OldLine = OldLine.AsLineEditable;
            if isempty(OldLine)
                error('Line %s is not editable',Line_Name);
            end
            OldLine.OverwriteWith(Line);
            Lines.Delete(Line);
        end
    end
end

function last_logged = print_new_log(EvApp, last_logged)
    % print out echoview log message if it has changed 
    logged = EvApp.GetLastLogMessage();
    if ~strcmp(logged, last_logged)
        last_logged = logged;
        fprintf('%s\n', logged);
    end
end
