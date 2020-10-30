function file_sets = generate_filelists(file_list, time_block)
% generate_filelists
%
% file_list is a list of Simrad raw file containing the date and time in
% the form -D20110920-T152230 as part of the file name.
%
% file_sets is a list of file lists, each list containing files which start
% within time_block hours of each other.
%
% Assumes file_list contains files in chronological order.
%
% The last file in each file list appears as the first file in the next
% file list.

%#ok<*AGROW>

if iscell(file_list)
    file_lists = file_list;
else
    file_lists{1} = file_list;
end

if exist(file_lists{1}, 'file') ~= 2
    warning('generate_filelists:file_not_found', 'file not found: %s', file_lists{1})
    file_sets{1,1} = {};
    return;
end

for i = 1:length(file_lists)
    if ~isempty(file_lists{i})
        % get the list of filenames
        fid = fopen(file_lists{i},'r');
        datafilelist = textscan(fid,'%q', ...
            'commentStyle', '#', ...
            'delimiter', '');
        fclose(fid);
        
        if time_block <= 0
            file_sets(1,i) = datafilelist;
            
        else
            % break into blocks of up to time_block hours
            
            datafilelist = datafilelist{1};
            
            % get start time of each file from file name
            file_time = NaN(size(datafilelist));
            for f = 1:length(datafilelist)
                try
                    timestamp = simrad_date_string(datafilelist{f});
                    timestamp = timestamp(regexp(timestamp, '[0-9]'));
                    
                    file_time(f) = datenum(timestamp, 'yyyymmddHHMMSS');
                catch e                                                         %#ok<NASGU>
                    warning('generate_filelists:filelist', 'Couldn''t include file %s in file list', datafilelist{f});
                end
            end
            
            % drop bad files
            datafilelist(isnan(file_time)) = [];
            file_time(end+1) = 99999999;   % add dummy end of the world
            
            if i == 1
                % break file_list up into blocks by time.
                file_set_id = 0;
                start = 1;
                hours = time_block / 24;
                
                for f = 2:length(file_time)
                    if file_time(f) - file_time(start) > hours
                        file_set_id = file_set_id + 1;
                        file_sets{file_set_id,i} = datafilelist(start:f-1);
                        start_time(file_set_id) = file_time(start);
                        finsh_time(file_set_id) = file_time(f);
                        start=f-1;
                    end
                end
            else
                % Find file set which covers each base set
                for f = 1:file_set_id
                    start = max(1, find(file_time >= start_time(f),1) - 1);
                    finsh = max(1, find(file_time >= finsh_time(f),1) - 1);
                    file_sets{f,i} = datafilelist(start:finsh);
                end
            end
        end
    end
end


