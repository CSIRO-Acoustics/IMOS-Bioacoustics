function UTC_adjust(indir, outdir, newnames, threshold, adjust, limit)
% UTC_ADJUST reads all the .raw files in indir and copies them to outdir
% with the record timestamp adjusted according to the UTC offset calculated
% from NMEA strings.
%
% Adjustments are calculated from GPS strings in the NME0 telegrams. No
% attempt is made to compensate for PC clock drift. Jumps in time exceeding
% the threshold level [default 5 minutes] cause the adjustment to be
% recalculated.
%
% INPUTS
%   indir       Directory containing raw files to process (all .raw files
%               in the directory will be processed
%               or a file containing a list of files to be be processed 
%               or a cell array of filenames to be processed
%               or if empty the user will be asked for the directory [].
%   outdir      Directory to write files to, if empty the user will be
%               asked to provide the directory [].
%   newnames    if true outfiles will be named according to the timestamp
%               of the first output record, otherwise infile names are used
%   theshold    Value used to detect a jump in time or need to adjust
%               times, offsets less than threshold will not be changed. 
%               [5 minutes]
%   adjust      Initial time adjustment. Adjustment to be made until a GPS
%               time is found in a NME0 telegram. If empty (default) data
%               will be buffered until a GPS record is read. This value
%               should only be supplied if the file lacks GPS time or there
%               is so much data before the first GPS time that buffering is
%               too resource intensive.
%   limit       Maximum adjustment, larger jumps are considered bad data.
%
% Author    Gordon.Keith@csiro.au
% Date      20141105


calendar0 = datenum(1601,1,1);      % Windows start of date numbers

if nargin < 1 || isempty(indir)
    indir = uigetdir('','Directory of .raw files to read');
end

if nargin < 2 || isempty(outdir)
    outdir = uigetdir(indir, 'Directory to write .raw files to');
elseif isempty(fileparts(outdir))
    outdir = fullfile(indir,outdir);
end

if nargin < 3 || isempty(newnames)
    newnames = false;
end

if nargin < 4 || isempty(threshold)
    threshold = 5/(24*60);           % 5 minutes
end

if nargin < 5
    adjust = [];
end

if nargin < 6
    limit = 1;
end

if ~exist(outdir, 'dir')
    mkdir(outdir)
end

if iscell(indir)
    filelist = indir;
elseif exist(indir, 'file') == 2
    fid = fopen(indir);
    scan = textscan(fid,'%[^\n]');
    filelist = scan{1};
else
    list = dir(fullfile(indir,'*.raw'));
    filelist = fullfile(indir,{list.name});
end

last = [];
lastadjust = 0;

buffer=struct();        % buffer data until we find a GPS time to calculate adjustment
buffer(1)=[];


for i = 1: length(filelist)
    infile = filelist{i};
    [~,filename,ext] = fileparts(filelist{i});
    outfile = fullfile(outdir,[filename,ext]);
    
    inf = java.io.File(infile);
    otf = java.io.File(outfile);
    if inf.isequal(otf)
        error('Infile and outfile are the same: %s!', infile);
    end
    
    fprintf('%s %s\n',datestr(now), filename);
    fid = fopen(infile);
    if newnames
        fout = 0;
    else
        fout = fopen(outfile,'w');
    end
    
    try
    while ~feof(fid)
        len = fread(fid,1,'uint32');
       
        if isempty(len)
            break;
        end
        
        if (len) > 1e6 || len < 12
            error('Bad record length')              % if matlab gets too large a read request it is uninterruptable
        end
        
        % read next telegram
        buffer(end+1).len = len;                    %#ok<AGROW>
        buffer(end).datagramtype = fread(fid,4,'*char')';
        
        buffer(end).lowdatetime  = fread(fid,1,'uint32');
        buffer(end).highdatetime = fread(fid,1,'uint32');
        buffer(end).datetime = (buffer(end).highdatetime*2^32 + buffer(end).lowdatetime) ...
            / 8.6400e+11 + calendar0;
        
        buffer(end).text = fread(fid,len - 12, '*char')';
        
        if abs(buffer(end).datetime - last) > threshold
            % jump in data, assume a new adjustment is needed.
            % apply last found adjustment to any buffered data.
            next = buffer(end);
            buffer(end) = [];
            fout = check_open(fout, buffer, outdir, lastadjust);
            writebuffer(fout, buffer, lastadjust);
            buffer = next;
            
            adjust = [];
        end
        last = buffer(end).datetime;
        
        if strcmp(buffer(end).datagramtype, 'NME0')
            nmea = parsenmea(buffer(end).text);
            if isfield(nmea,'time')
                offset = nmea.time - buffer(end).datetime;
                if nmea.time < 1    % does not includes date
                    offset = mod(offset,1);
                    if offset > 0.6
                        offset = offset -1;
                    end
                end
                
                if abs(offset) > limit
                    warning('At %s offset %s greater than %g days', ...
                        datestr(buffer(end).datetime),timestr(offset), limit);
                end
                
                if isempty(adjust)
                    adjust = 0;
                end
                
                if abs(adjust - offset) > threshold && abs(offset) < limit
                    if adjust ~= 0
                        warning('Change in offset from %s to %s', ...
                            timestr(adjust), timestr(offset))
                    end
                    adjust = offset;
                end
                
                lastadjust = adjust;
            end
        end
        
        if ~isempty(adjust)
            % if we have a current adjustment, flush the buffer.
            fout = check_open(fout, buffer, outdir, adjust);
            buffer = writebuffer(fout, buffer, adjust);
        end
        
        if len ~= fread(fid,1,'uint32')
            error('Corrupt file, record length mismatch');
        end

    end
    catch e
        warning('Problem reading raw file \n%s\n%s', infile, e.message);
    end
    % flush buffer and close file
    fout = check_open(fout, buffer, outdir, lastadjust);
    writebuffer(fout, buffer, lastadjust);
    fclose(fid);
    fclose(fout);
end

    function buffer=writebuffer(fout, buffer, adjust)
    % write all records in the buffer to fout applying time adjustment
    % adjust.
        while ~isempty(buffer)
            if adjust ~= 0
                buffer(1).datetime = buffer(1).datetime + adjust;
                dt = uint64((buffer(1).datetime - calendar0) * 8.6400e+11);
                buffer(1).lowdatetime = uint32(bitand(dt,hex2dec('0ffffffff')));
                buffer(1).highdatetime = uint32(bitand(bitshift(dt,-32),hex2dec('0ffffffff')));
            end
            
            fwrite(fout,buffer(1).len,'uint32');
            fwrite(fout,buffer(1).datagramtype,'char');
            fwrite(fout,buffer(1).lowdatetime,'uint32');
            fwrite(fout,buffer(1).highdatetime,'uint32');
            fwrite(fout,buffer(1).text,'char');
            fwrite(fout,buffer(1).len,'uint32');

            buffer(1) = [];
        end
    end

    function fout = check_open(fout, buffer, outdir, adjust)
    % if fout is not a valid file descriptor open a new file with name
    % based on the timestamp of the first record in buffer.
    
        if fout > 0 || isempty(buffer)
            return;
        end
        
        outf = sprintf('D%s.raw',datestr(buffer(1).datetime + adjust, 'yyyymmdd-THHMMSS'));       
        fout = fopen(fullfile(outdir,outf),'w');
        fprintf('%s-%s\n', datestr(now), outf);
    end

    function str = timestr(time)
    % A string given signed time
    
        tim = sprintf('%+f', time);
        dot = find(tim == '.', 1, 'first');
        tim(dot) = ' ';
        if tim(2) == '0'
            dot=1;
        end
        tim(dot+1:end) = [];
        str = sprintf('%s%s', tim, datestr(abs(time),'HH:MM:SS'));
    end
end


