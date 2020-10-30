function [listFile,infFile,gpsFile,pitchFile,rollFile,depthFile] = getGPS(infiles, outfile, varargin)
% getGPS takes one or more EK60 compatible files and generates:
% outfile_data_files.txt - list of files if length(infiles > 1)
% outfile.inf -
% outfile.gps.csv
% outfile.pitch.csv
% outfile.roll.csv
% outfile.depth.csv
%
% Inputs
%       infiles - cell array of file names
%       outfile - root name of output files
%

% get list of files if not provided
if nargin < 1
    infiles='';
end

if isempty(infiles) || ~iscell(infiles)
       
    [files,path] = uigetfile('*.raw','Select EK60 files',infiles,'MultiSelect','on');
    if iscell(files)
        infiles = fullfile(path,sort(files));
    elseif files == 0
        error('No files specified');
    else
        infiles = cell(1);
        infiles{1} = fullfile(path,files);
    end
end

% get outfile if not provided
if nargin < 2 
    outfile = '';
end
if isempty(outfile) || exist(outfile,'dir') == 7
    [file,path] = uiputfile('*','Specify the base name for output files',outfile);
    outfile = fullfile(path,file);
end

% check for readdatagram
if isempty(which('readdatagram'))
    root_path = fileparts(mfilename('fullpath'));
    if exist(fullfile(root_path, 'read_write_ek60'),'dir') == 7
        addpath(fullfile(root_path, 'read_write_ek60'));
    else
        root_path = uigetdir(root_path,'Directory containing readdatagram.m');
        if ~isequal(root_path,0)
            addpath(root_path);
        end
    end
end

% check arguments
% info = true;        % if true generate .inf file otherwise don't generate stats
info = false;        % if true generate .inf file otherwise don't generate stats
                    % set info to false as code was crashing. Tim Ryan 2018
                    % 03 
existing = false;   % use existing .inf and .gps.csv files if they exist
detail = 'rawheader';

for i = 3:nargin
    if strncmpi(varargin{i-2}, 'noinf', 5)
        % don't generate .inf file, around 5% faster
        info = 0;
        detail = 'noraw';
    end
    if strncmpi(varargin{i-2}, 'existing', 3)
        existing = true;
    end
end

% output files, 
listFile  = '';
infFile   = '';
gpsFile   = '';
pitchFile = '';
rollFile  = '';
depthFile = '';
list  = -1;
inf   = -1;
gps   = -1;
pitch = -1;
roll  = -1;
depth = -1;

if length(infiles) > 1
    listFile = [outfile '_data_files.txt'];
    list  = fopen(listFile,'w');
    if list < 1 
        error(['Unable to open list file:',listFile]);
    end
end

if info
    infFile = [outfile '.inf'];
    inf   = fopen(infFile,'w');
end

% statistics variables, used to write .inf file
gpscount=0;     % number of GPS fixes
count = [];     % number of pings per channel
settings = {};  % channel power and pulse changes
minlat = 999;   % southern most latitude
maxlat = -999;  % northern most latitude
minlon = 999;   % western most longitude
maxlon = -999;  % eastern most longitude
start = [];     % first GPS fix
last = [];      % last GPS fix

offsets = containers.Map(0,0);
offset0 = [];
jitter = 1;

% process files
%<<<<<<< .mine
%for f = 1:length(infiles)    
%    keyboard
% ||||||| .r845
%for f = 1:length(infiles)
%    in = fopen(infiles{f});

% Added these lines as infiles was coming in as a single 1 b 1 cell. Not
% sure how it ever work. Commented out 09/03/2018 as unsure if there will
% be unintended consequences with other proceessing. 
% % % % %     [m,n] = size(infiles)
% % % % %     if isequal(m,1) & isequal(n,1)
% % % % %         for i=1:length(infiles{1})
% % % % %             newinfiles{i} = infiles{1}{i};
% % % % %         end
% % % % %     end
% % % % %     infiles = newinfiles;


%=======
for f = 1:length(infiles)
%>>>>>>> .r1033

    fprintf('%s %d/%d: %s\n', datestr(now,'HH:MM:SS'), f, length(infiles), infiles{f});    
    if list > 0
        fprintf(list, '%s\n', infiles{f});
    end
    
    if existing && exist([infiles{f} '.gps.csv'], 'file') == 2
        if info
            warning('GPS:inf_exist', 'Use of existing .inf files not yet supported')
        end
        
        if gps == -1
            gpsFile = [outfile '.gps.csv'];
            gps   = fopen(gpsFile,'w');
            fprintf(gps,'GPS_date,GPS_time,GPS_milliseconds,Latitude,Longitude,UTC_offset\n');
        end
        gpsin = fopen([infiles{f} '.gps.csv']);
        fgetl(gpsin);
        while ~feof(gpsin)
            fprintf(gps,'%s\n',fgetl(gpsin));
        end
        fclose(gpsin);        
        if exist([infiles{f} '.pitch.csv'], 'file') == 2
            if pitch < 0
                % open the files the first time we want to use them
                pitchFile = [outfile '.pitch.csv'];
                rollFile = [outfile '.roll.csv'];
                pitch = fopen(pitchFile,'w');
                roll  = fopen(rollFile,'w');
                fprintf(pitch,'Pitch_date,Pitch_time,Pitch_milliseconds,Pitch_angle\n');
                fprintf(roll,'Roll_date,Roll_time,Roll_milliseconds,Roll_angle\n');
            end            
            pitchin = fopen([infiles{f} '.pitch.csv']);
            fgetl(pitchin);
            while ~feof(pitchin)
                fprintf(pitch,'%s\n',fgetl(pitchin));
            end
            fclose(pitchin);
            
            rollin = fopen([infiles{f} '.roll.csv']);
            fgetl(rollin);
            while ~feof(rollin)
                fprintf(roll,'%s\n',fgetl(rollin));
            end
            fclose(rollin);
        end
        if exist([infiles{f} '.depth.csv'], 'file') == 2
            if depth < 0
                depthFile = [outfile '.depth.csv'];
                depth = fopen(depthFile,'w');
                fprintf(depth,'Depth_date,Depth_time,Depth_milliseconds,Platform_Depth\n');
            end
            
            depthin = fopen([infiles{f} '.depth.csv']);
            fgetl(depthin);
            while ~feof(depthin)
                fprintf(depth,'%s\n',fgetl(depthin));
            end
            fclose(depthin);
        end
    else
    in = fopen(infiles{f});
    
    datagram = readdatagram(in);
    
    if ~strcmp(datagram.datagramtype,'CON0')
        error(['Expecting config telegram at start of ' infiles{f}]);
    end
    
    con = datagram.CON0;
    if con.transducercount ~= length(count)
        if ~isempty(count)
            error(['Incompatible files - transducer count changes \n' ...
                num2str(length(count)) ': ' infiles{f-1} '\n' ...
                num2str(con.transducercount) ': ' infiles{f} ]);
        end
        
        count = zeros(con.transducercount,1);
        power = zeros(con.transducercount,1);
        pulse = zeros(con.transducercount,1);
    end
    
    datagram = readdatagram(in,detail);
    while ~isempty(datagram)
        if info && strcmp(datagram.datagramtype, 'RAW0')
            % get settings and counts for .inf file.
            channel = datagram.RAW0.channel;
            count(channel) = count(channel) +1; %#ok<AGROW> it doesn't the analyser is confused
            
            if power(channel) ~= datagram.RAW0.transmitpower || ...
               pulse(channel) ~= datagram.RAW0.pulselength
                settings{end+1} = ...
                    ['[' trim(con.transducer(channel).channelid) ']  (' ...
                    datestr(datagram.datetime, 'yyyy-mm-dd HH:MM:SS') ')  ' ...
                    num2str(datagram.RAW0.pulselength * 1000000) ' us   ' ...
                    num2str(datagram.RAW0.transmitpower) ' W']; %#ok<AGROW>
                power(channel) = datagram.RAW0.transmitpower;
                pulse(channel) = datagram.RAW0.pulselength;
            end
        end
        
        if strcmp(datagram.datagramtype, 'NME0')
            if isfield(datagram.NME0, 'latitude')
                
                % calculate .inf file statistics
                if info
                    gpscount = gpscount + 1;
                    if isempty(start)
                        start = datagram;
                        distance = 0;
                        last = datagram;
                    end
                    if datagram.datetime - last.datetime < 1/48     %don't count gaps bigger than half an hour
                        distance = distance + dist(last.NME0, datagram.NME0);
                    end
                    last = datagram;
                    
                    if minlat > datagram.NME0.latitude ; minlat = datagram.NME0.latitude ; end
                    if maxlat < datagram.NME0.latitude ; maxlat = datagram.NME0.latitude ; end
                    if minlon > datagram.NME0.longitude ; minlon = datagram.NME0.longitude ; end
                    if maxlon < datagram.NME0.longitude ; maxlon = datagram.NME0.longitude ; end
                end
                
                % calculat UTC offset if GPS has time
                offset = [];
                if isfield(datagram.NME0, 'time')
                    offset = datagram.datetime - datagram.NME0.time;
                    offset = offset - round(offset);
                    offset = offset * 24 * 60 * 60;
                    
                    if isempty(offset0)
                        offset0 = offset;
                    end
                    off = round((offset - offset0)/jitter);
                    if offsets.isKey(off)
                        offsets(off) = offsets(off) + 1;
                    else
                        offsets(off) = 1;
                    end
                end
                
                if gps == -1
                    gpsFile = [outfile '.gps.csv'];
                    gps   = fopen(gpsFile,'w');
                    fprintf(gps,'GPS_date,GPS_time,GPS_milliseconds,Latitude,Longitude,UTC_offset\n');
                end
                
                printtime(gps,datagram.datetime);
                fprintf(gps,'%g,%g', datagram.NME0.latitude, datagram.NME0.longitude);
                if ~isempty(offset)
                    fprintf(gps,',%.3f', offset);
                end
                fprintf(gps,'\n');
            end
            
            if strcmp(datagram.NME0.sentence, 'PRY')
                % write pitch and roll .csv files
                
                if pitch < 0
                    % open the files the first time we want to use them
                    pitchFile = [outfile '.pitch.csv'];
                    rollFile = [outfile '.roll.csv'];
                    pitch = fopen(pitchFile,'w');
                    roll  = fopen(rollFile,'w');
                    fprintf(pitch,'Pitch_date,Pitch_time,Pitch_milliseconds,Pitch_angle\n');
                    fprintf(roll,'Roll_date,Roll_time,Roll_milliseconds,Roll_angle\n');
                end
                
                printtime(pitch,datagram.time);
                fprintf(pitch,'%.2f\n',datagram.NME0.pitch);
                printtime(roll,datagram.time);
                fprintf(roll,'%.2f\n',datagram.NME0.roll);
                
            end
            
            if isfield(datagram.NME0,'depth')
                if depth < 0
                    depthFile = [outfile '.depth.csv'];
                    depth = fopen(depthFile,'w');
                    fprintf(depth,'Depth_date,Depth_time,Depth_milliseconds,Platform_Depth\n');
                end
                
                printtime(depth,datagram.time);
                fprintf(depth,'%.2f\n',datagram.NME0.depth);              
                
            end
        end
        
        datagram = readdatagram(in,detail);
    end
    
    fclose(in);
    end
end

if info
    % write .inf file
    fprintf(inf,'\nTrack Data File:     ');
    [~, name, ext] = fileparts(outfile);
    fprintf(inf, [name ext]);
    fprintf(inf,'\nTrack Data Path:     ');
    fprintf(inf, '%s\n',outfile);
    keyboard    
    fprintf(inf,'\n\nChannels:     \n  ');
    for i = 1:con.transducercount
        fprintf(inf,trim(con.transducer(i).channelid));
        fprintf(inf,' : %d\n  ', count(i));
    end
    
    fprintf(inf,'GPS fixes : %d\n', gpscount);
    
    if ~isempty(offset0)
        if offsets.Count == 1
            fprintf(inf,'    UTC offset : %d seconds\n', round(offset0/jitter) * jitter);
        else
            fprintf(inf,'    UTC offsets :\n');
            off = offsets.keys;
            for i = 1:length(off)
                fprintf(inf,'       %d seconds : %d\n', round((offset0+off{i})/jitter) * jitter, offsets(off{i}));
            end
            if off{end} - off{1} > 300
                warning('GPS:OFFSET','UTC offset varies from %d to %d seconds in %s', ...
                    off{1},off{end}, outfile);
                warndlg(sprintf('UTC offset varies from %d to %d seconds in \n%s', ...
                    off{1},off{end}, outfile), 'Large time jump');
            end
        end
    end
    
    for i =1:length(settings);
        fprintf(inf,'\n        ');
        fprintf(inf,settings{i});
    end
    
    duration = (last.datetime - start.datetime) * 24;
    speed = distance / duration;
    fprintf(inf,'\n\nNavigation Totals:\nTotal Time:         % 9.4f hours\n', duration);
    fprintf(inf,'Total Track Length: % 9.4f km\n', distance *1.852);
    fprintf(inf,'Average Speed:    % 9.4f km/hr (%.4f knots)\n', speed * 1.852, speed);
    
    fprintf(inf,'\nStart of Data:\nTime:  ');
    fprintf(inf,datestr(start.datetime,'yyyy-mm-dd HH:MM:SS.FFF'));
    dv=datevec(start.datetime);
    dv(2:6)=0;
    fprintf(inf,'  JD%d\n',floor(start.datetime - datenum(dv)));
    fprintf(inf,'Lon:  % 9.4f    Lat: % 9.4f\n', start.NME0.longitude, start.NME0.latitude);
    
    fprintf(inf,'\nEnd of Data:\nTime:  ');
    fprintf(inf,datestr(last.datetime,'yyyy-mm-dd HH:MM:SS.FFF'));
    dv=datevec(last.datetime);
    dv(2:6)=0;
    fprintf(inf,'  JD%d\n',floor(last.datetime - datenum(dv)));
    fprintf(inf,'Lon:  % 9.4f    Lat: % 9.4f\n', last.NME0.longitude, last.NME0.latitude);
    
    fprintf(inf,'\nLimits\n');
    fprintf(inf,'Minimum Longitude:     % 9.4f   Maximum Longitude:     % 9.4f\n', minlon, maxlon);
    fprintf(inf,'Minimum Latitude:      % 9.4f   Maximum Latitude:      % 9.4f\n', minlat, maxlat);
    fclose(inf);
end

if gps > 0
    fclose(gps);
end 
if list > 0
    fclose(list);
end
if pitch > 0
    fclose(pitch);
    fclose(roll);
end
if depth > 0
    fclose(depth);
end



function printtime(fid, time)
% Print time stamp in echoview csv format
fprintf(fid,datestr(time, 'yyyy-mm-dd,HH:MM:SS,FFF,'));

function distance = dist(from, to)
% returns approximate distance in Nm
% upgrade to a more accurate calculation if you want better results

distance = sqrt((to.latitude - from.latitude) ^2 + ...
    (cos(to.latitude * pi() / 180) * (to.longitude - from.longitude)^2)) * 60;


function text=trim(text)
text(text == 0) = [];
text=strtrim(text);
