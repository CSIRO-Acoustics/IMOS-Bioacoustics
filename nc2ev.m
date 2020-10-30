function nc2ev(ncfile)
% function to convert variables named in the matrix_variable_list to Echoview
% .sv.csv format
% stage 1: export names variable to echoview ascii format
% stage 2: import into an Ev file
%
% Preconditions:
%   full path to variable ncfile is given
%   directory for storing the Echoview outputs will be created based on the
%   path to the variable ncfile. Evfile will also be created in that
%   directory
%   call viz_sv to extract the data from the netCDF file
%
% useage:
% nc2ev(ncfile) where ncfile is the full path to the netCDF file.
% Tim Ryan 31/07/2018

%% Extract data from ncfile using the utility ncfile (note this function may eventually be built into viz_sv in which case this step would be removed.

    if isequal(nargin,1)
        % default exports
        matrix_variable_list = [{'Sv'}    {'Svraw'}    {'signal_noise'}    {'pg'}]; % i.e. echograms and the like
        vector_variable_list = [{'NASC'},{'NASCraw'},{'background_noise'},{'epipelagic'},{'upper_mesopelagic'},{'lower_mesopelagic'},{'day'}]; % i.e. lines and the like, background noise, epipelagic etc 
    end    
    template_dir = 'Q:\IMOS_BASOOP\Echoview_NetCDF'; % set this as location to find echoview templates for netCDF format data. 

%% Create a folder for the echoview ascii outputs, export if needed, create ev file from exports

    [dir_root, file] = fileparts(ncfile); % changing to 'dir_root' as it conflict with matlab inbuilt 'dir'
    if ~isdir([dir_root '\echoview_ascii'])
        mkdir([dir_root '\echoview_ascii']);
    end
    cd([dir_root '\echoview_ascii'])
    
    file_name = ncreadatt(ncfile, '/', 'deployment_id'); % reading global attribute 'deployment_id' for csv, ev file naming convention
    
    % Haris 11 July 2019 - I tried to create file names based on the NetCDF
    % file name. However, due to some reason Echoview is not allowing to
    % add those into the filesets. May be due to length of file name. Now
    % using unique 'deployment_id' to prefix file names of csv, evl, evd,
    % and ev file. This would be useful if we want to submit these files to
    % AODN and unique file name is needed for them.

    csv_info = dir('*.csv'); % making a list of CSV files to check 
    ev_info = dir('*.ev');   % making a list of ev files to check 
    
    if isempty(csv_info)
        [data] =  export_to_ascii(ncfile,dir_root,matrix_variable_list,vector_variable_list);
    else 
        redo = questdlg('Echoview ASCII files already exist, do you want to re-export?','Redo export','Yes','No','No'); % default to 'No' because it already exist
        if isequal(redo,'Yes')
            delete * 
            % Haris 11 July 2019 - delete all files in the folder assuming
            % that for redo export all files need to be new (csv, evl, evd,
            % and ev). Note files have new file name prefix with
            % deployment_id.
            [data] = export_to_ascii(ncfile,dir_root,matrix_variable_list,vector_variable_list); % call function to export ASCII files and evl 
        end
    end 
    
    if ~exist('data') % i.e. did not generate this in the previous stepm then will need to rat it out using viz_sv;
        [data] = viz_sv(ncfile,[],'noplots','all');
    end
    
    if isempty(ev_info)
        create_ev_from_ascii([dir_root '\echoview_ascii'],matrix_variable_list, vector_variable_list,template_dir,data); % call function to create ev files from exports 
    else
        redoev = questdlg('Echoview file already exist, do you want to re-create?','Redo ev file','Yes','No','No'); 
        if isequal(redoev,'Yes')                        
            delete *.ev *.evd  
            % Haris 11 July 2019 - delete all Echoview related stuff -
            % assumption is for redo ev new ev file will be created
            % prefixing global attribute 'deployment_id'
            create_ev_from_ascii([dir_root '\echoview_ascii'],matrix_variable_list, vector_variable_list,template_dir,data);     
        end
    end    

%% Function to export ASCII files from the NetCDF

    function [data] = export_to_ascii(ncfile,dir_root,matrix_variable_list,vector_variable_list)
        
        [data] = viz_sv(ncfile,[],'noplots','all');       

        ns = length(data.depth);
        h = (data.depth(end) - data.depth(1)) / (ns-1);
        range = sprintf('%g, %g', data.depth(1) - h/2, data.depth(end) + h/2);
        
        fprintf('                     exporting ASCII files to %s\n',[dir_root '\echoview_ascii']) % to keep aligned with process_BASOOP progress messages
        
        for j=1:length(matrix_variable_list)    
            for k=1:length(data.channels)        
                if isequal(length(data.channels),1) % single freq            
                    if eval(['min(size(data.' matrix_variable_list{j} '))']) > 1                        
                       output_file = [matrix_variable_list{j} '_' data.channels{k} '.sv.csv'];
                       output_file = strcat(file_name,'_',output_file); % prefix global attribute 'deployment_id'
                       fid = fopen(output_file,'w+');
                       output_var = matrix_variable_list{j};
                    end            
                elseif eval(['length(size(data.' matrix_variable_list{j}, '))'])>2                                 
                       output_file = [matrix_variable_list{j} '_' data.channels{k} '.sv.csv'];
                       output_file = strcat(file_name,'_',output_file); % prefix global attribute 'deployment_id'
                       fid = fopen(output_file,'w+');
                       output_var = matrix_variable_list{j};
                end                                     
                fprintf(fid,'Ping_date,Ping_time,Range_start,Range_stop,Sample_count\n');
                for i = 1:length(data.latitude)            
                    dt = datestr(data.time(i), 'yyyy-mm-dd, HH:MM:SS');           
                    Output_data = eval(['{data.' output_var '(:,i,k)};']);                                                         
                    Output_dte = sprintf('%.1f,',Output_data{:});
                    Output_dte = strrep(Output_dte,'-Inf','-999');                    
                    fprintf(fid,'%s,%s,%d,',...
                        dt,range,ns);            
                    fprintf(fid,'%s',Output_dte(1:end-1));
                    fprintf(fid,'\n');               
                end
%                 fprintf('exported variable %s\n to %s''\''%s\n',dir_root,output_var, output_file)
                fprintf('                     exported variable %s\n',output_var)
            end
        end
        
        % now print out a gps.csv file to same location
        
        output_file_gps = strcat(file_name,'_','position.gps.csv'); % prefix global attribute 'deployment_id'
        gid = fopen(output_file_gps,'w+');        
%         gid = fopen([dir_root '\echoview_ascii\position.gps.csv'],'w+');
        fprintf(gid, 'GPS_date, GPS_time, Latitude, Longitude\n');
        for i = 1:length(data.latitude)
            dt = datestr(data.time(i), 'yyyy-mm-dd, HH:MM:SS');
            fprintf(gid, '%s, %g, %g\n', dt, data.latitude(i), data.longitude(i));
        end
        
        fprintf('                     exported variable position\n') % to keep aligned with process_BASOOP progress messages        
        fprintf('                     starting line export\n'); 

        for k = 1:length(vector_variable_list)
            if (length(data.channels)>1)
                for i=1:length(data.channels)
                    vec_var = ['data.' vector_variable_list{k}];
                    if isequal(eval(['length(size(data.' vector_variable_list{k} '))']),3) % we have the three dimensional matrix
%                         eval(['fid = fopen(' '''' pwd '\' vector_variable_list{k} '_' data.channels{i}  '.evl'', ''w+'');'])
                        eval(['fid = fopen(' '''' pwd '\' file_name '_' vector_variable_list{k} '_' data.channels{i}  '.evl'', ''w+'');']) % prefix global attribute 'deployment_id'
                        datastr = [vec_var '(1,j,i)'];
                    elseif length(eval(['data.' vector_variable_list{6} '(:,1)' ])) > 1  & ~isequal(vector_variable_list{k}, 'day')
%                         eval(['fid = fopen(' '''' pwd '\' vector_variable_list{k} '_' data.channels{i}  '.evl'', ''w+'');'])
                        eval(['fid = fopen(' '''' pwd '\' file_name '_' vector_variable_list{k} '_' data.channels{i}  '.evl'', ''w+'');']) % prefix global attribute 'deployment_id'
                        datastr = [vec_var '(i,j)'];
                    else  % one dimensional matrix (e.g. day)
%                         eval(['fid = fopen(' '''' pwd '\' vector_variable_list{k} '.evl'', ''w+'');'])
                        eval(['fid = fopen(' '''' pwd '\' file_name '_' vector_variable_list{k} '.evl'', ''w+'');']) % prefix global attribute 'deployment_id'
                        datastr = [vec_var '(j)'];
                    end
                    if ~isequal(vector_variable_list{k},'day')
                        s = size(eval(vec_var));
                        len_dataset = s(2);
                        n = len_dataset;
                    else
                        s = size(eval(vec_var));
                        len_dataset = s(1);
                        n = len_dataset;
                    end
                    % this loop tallys number of valid rows as there as some inf's
                    % kicking around
                    for j=1:n
                        value = eval(datastr);
                        if isinf(value) | isnan(value)
                            len_dataset = len_dataset-1;
                        end
                    end
                    fprintf(fid,'EVBD 3 3.00.41\n%d\n',len_dataset);
                    for j=1:n % i.e. this is the length of the data
                        dt = datestr(data.time(j), 'yyyymmdd HHMMSS0000');
                        value = eval(datastr);
                        if (~isinf(value) | ~isnan(value)) & ~isequal(vector_variable_list{k},'day') % there are some infs around - don't print that line
                            fprintf(fid, '%s %f 1\n',dt, value);
                        end
                        if (~isinf(value) | ~isnan(value)) & isequal(vector_variable_list{k},'day') % there are some infs around - don't print that line
                            fprintf(fid, '%s %f 1\n',dt, value*100);
                        end
                    end
                    fclose(fid);
                end
            else % handle instance for single frequency
                vec_var = ['data.' vector_variable_list{k}];
                datastr = [vec_var '(j)'];
%                 eval(['fid = fopen(' '''' pwd '\' vector_variable_list{k} '.evl'', ''w+'');'])
                eval(['fid = fopen(' '''' pwd '\' file_name '_' vector_variable_list{k} '.evl'', ''w+'');']) % prefix global attribute 'deployment_id'
                n = length(data.NASC);
                len_dataset =n;
                for j=1:n
                    value = eval(datastr);
                    if isinf(value) | isnan(value)
                        len_dataset = len_dataset-1;
                    end
                end
                fprintf(fid,'EVBD 3 3.00.41\n%d\n',len_dataset);
                for j=1:n % i.e. this is the length of the data
                    dt = datestr(data.time(j), 'yyyymmdd HHMMSS0000');
                    value = eval(datastr);
                    if ~isinf(value) & ~isnan(value)
                        fprintf(fid, '%s %f 1\n',dt, value);
                    end
                end
            end
        end
        fprintf('                     finished line exports\n');
        fprintf('                     finished exporting NetCDF contents to Echoview ASCII format\n');        
        fclose all;
    end    
    
%% function to create ev files from ASCII exports

    function create_ev_from_ascii(asciidir,matrix_variable_list,vector_variable_list, template_dir,data)
        % use com scripting to create an ev file    
        S = dir(fullfile(asciidir,'*.csv'));      

        csvfiles = natsortfiles({S.name});
        %csvfiles = ls([asciidir '\*.csv']);    
        EvApp = actxserver('EchoviewCom.EvApplication');    
        EvFile = EvApp.NewFile;                    
        % sort the csv files for freq and place position variable at the end. 
        new_csvlist = []; position_files = [];    
        for i=1:length(csvfiles)
            csvfile = csvfiles(i);
            if strfind(cell2mat(csvfile),'position')
                position_file = cell2mat(csvfile);
            else
                new_csvlist = [new_csvlist;csvfile];
            end        
        end
        csvfiles = [new_csvlist; position_file];        
        for i=2:length(csvfiles)
            EvFile.Filesets.Add;        
        end
        for i=1:length(csvfiles)
            fileset2add = cell2mat(csvfiles(i));
            dt = strfind(fileset2add,'.');
            fileset2add = fileset2add(1:dt-1);
            EvFile.Filesets.Item(i-1).Name = fileset2add;
            datafile = [asciidir '\' cell2mat(csvfiles(i))];
            EvFile.Filesets.Item(i-1).DataFiles.Add(datafile);
        end    
        cd(asciidir); 
%         clipboard('copy', [asciidir '\ncdata.evd']) 
        evd_file = strcat(asciidir,'\',file_name,'.','evd'); % prefix global attribute 'deployment_id'
        clipboard('copy', evd_file)
        fprintf('                     go to the message box, follow instructions and then click OK button\n'); % to keep aligned with process_BASOOP progress messages 
        
        % Haris 12 July 2019 - to make message box and text bigger
                   
        instructions = ['\n=====================================================\n\n',...
            'Do the following steps in sequence:\n\n',...
            '(1)  Switch back to the Echoview window opened.\n\n',...
            '(2)  Open an echogram from the "Dataflow" window (does not matter which one).\n\n',...
            '(3)  From the top menu select: "Echogram > Export > To Echoview Data File Format".\n\n',...
            '(4)  When the dialog box opens "Select All" and then "Export". Do not change anything.\n\n',...
            '(5)  Export file name (*.evd) has been programmatically copied to the clipboard for you. Do not worry.\n\n',...
            '(6)  In the next dialog box, use keyboard shortcut "Control-V" for pasting clipboard copied "File name" and click "Save".\n\n',...
            '(7)  When saving is done, minimise Echoview window and click this "OK" button below.\n\n',...
            '====================================================\n'];
        
        h = msgbox(sprintf(instructions),'Help Dialog');
        set(h, 'position', [100 100 480 480]); % makes box bigger
        ah = get(h, 'CurrentAxes');
        ch = get(ah, 'Children');
        set(ch, 'FontSize', 13);  % makes text bigger
        uiwait(h)
        
        EvFile.Close; %exit this file

        template_sub_dir = 'F';
        for i=1:length(data.channels)
            template_sub_dir = [template_sub_dir '_' data.channels{i}];
        end
        template_sub_dir = strrep(template_sub_dir, 'kHz','');    
        EvFile = EvApp.OpenFile([template_dir '\' template_sub_dir '\nc_template.ev']);
        if isempty(EvFile)
             msgbox('Template file to create EV file does not exist - check ''Q:\IMOS_BASOOP\Echoview_NetCDF''','Template not found','warn');
             return
        else
            EvFile.Filesets.Item(0).DataFiles.RemoveAll;    % get rid to template dummy file
%             EvFile.Filesets.Item(0).Datafiles.Add([pwd '\ncdata.evd']); % add the previously exported evd file
            EvFile.Filesets.Item(0).Datafiles.Add(evd_file); % evd file name now include global attribute 'deployment_id' 
        end    
        % add the evlfiles according to named vector variable list. 
        fprintf('                     importing line data to the EV file\n')
        
        if (length(data.channels)>1)
            for k = 1:length(vector_variable_list)
                if isequal(vector_variable_list{k},'day')
                    try
%                         EvFile.Import([pwd '\' vector_variable_list{k} '.evl']);
                        EvFile.Import([pwd '\' file_name '_' vector_variable_list{k} '.evl']); % prefix global attribute 'deployment_id'
                    catch
                        fprintf('failed to import the line %s\n',[file_name '_' vector_variable_list{k} '.evl'])
                    end        
                else
                    for j=1:length(data.channels)
                        try
%                             EvFile.Import([pwd '\' vector_variable_list{k} '_' data.channels{j} '.evl']);
                            EvFile.Import([pwd '\' file_name '_' vector_variable_list{k} '_' data.channels{j} '.evl']); % prefix global attribute 'deployment_id'
                        catch
                            fprintf('failed to import the line %s\n',[file_name '_' vector_variable_list{k} '_' data.channels{j} '.evl'])
                        end
                    end
                end
            end    
        else
            for k = 1:length(vector_variable_list)
                try
%                     EvFile.Import([pwd '\' vector_variable_list{k} '.evl']);
                    EvFile.Import([pwd '\' file_name '_' vector_variable_list{k} '.evl']); % prefix global attribute 'deployment_id'
                catch
                    fprintf('failed to import the line %s\n',[file_name '_' vector_variable_list{k} '.evl'])
                end     
            end
        end
        
        nc_viewer = strcat(asciidir,'\',file_name,'.','ev'); % prefix global attribute 'deployment_id'
        
        EvFile.SaveAs(nc_viewer); 
        EvFile.Close;
        % Implemeted in the 'viz_sv' 
% %         msgbox(sprintf('Finished creating ev file "nc_viewer.ev"\nIt Echoview will now open file. Close when done\n')) 
% %         ques = questdlg('Do you want to open EV file for viewing NetCDF?','View NetCDF in Echoview','Yes','No','Yes');
% %         if isequal(ques,'Yes')
% %             winopen(nc_viewer)
%         end           
    end
end