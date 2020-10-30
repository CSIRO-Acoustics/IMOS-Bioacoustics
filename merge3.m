function filename = merge3( directory, outfile, frequency )
%MERGE3 merges *_Final_38kHz_cleaned.csv, *_HAC Sv 38 kHz.csv and
% _Reject_38kHz.csv into a single .csv file suitable for importing via the 
% echoview Parser in the IMOS-toolbox.
%
% Inputs:
%   directory - directory containing the input .csv files, the user will be
%   asked to provide a value if it not provided or empty.
%   outfile - full path the the output .csv file, the user will be asked to
%   provide a value if it is not provided or empty.
%   frequency - frequency of data to process this is part of the filename
%   of the input .csv files, if not provided 38 is used.
%
% Outputs:
%   filename - the full path to the output .csv file, this is the same as
%   outfile if provided.
%



% Depths below this will not be included in the output file.
MAX_DEPTH=1200;

% Pct_good below this will not be included in the output file.
MIN_GOOD=0;

% ensure the IMOS toolbox is in the Path

if isempty(strfind(path,pwd))
    addpath(genpath(pwd));

    % Add the ddb.jar java library and jdbc drivers to the classpath
    jars = dir(['Java' filesep '*.jar']);
    for j = 1 : length(jars)
        javaaddpath(['Java' filesep jars(j).name]);
    end
end

% Ask the user for the directory to use

if (nargin < 1 || isempty(directory))
    directory = uigetdir('q:\Processed_data');
end

% Ask user for the name of the output file

if (nargin < 2 || isempty(outfile))
    [file outpath] = uiputfile(fullfile(directory, 'merged.csv'));
    outfile = fullfile(outpath, file);
end

if (nargin < 3)
    frequency = 38;
end

channel = [ num2str(frequency) 'kHz'];

% Create list of *_Final_38kHz_cleaned.csv files

pattern = [ '^(.+)_Final_' channel '_cleaned\.csv$'];

sections = listFiles(directory, pattern);       % IMOS-toolbox/Util/listFiles.m

filename = outfile;
out = fopen(filename, 'w');

if out == -1
    error(['Can''t open ' filename]);
end

% read each file in order

first = 1;
pending = -99999;
currlayer = 0;
buffer = cell([1 100]);
intrval = 2;
layer = 3;
layer_min = 0;
layer_max = 0;
samples = 0;
latitude = 0;

sv_mean = 0;
sv_sd = 0;
sv_skew = 0;
sv_kurt = 0;

linelast = [];

for f = 1:length(sections)

  prefix = sections{f};

  cleanfile = fullfile(directory, [ prefix '_Final_' channel '_cleaned.csv' ]);
  rawfile = fullfile(directory,[ prefix '_HAC Sv ' channel '.csv' ]);
  countfile = fullfile(directory,[ prefix '_Reject_' channel '.csv' ]);

  if exist(rawfile,'file') ~= 2
      rawfile = fullfile(directory, [prefix '_HAC Sv 38 kHz.csv']);
  end
  
  clean = fopen(cleanfile, 'rt');
  raw = fopen(rawfile, 'rt');
  count = fopen(countfile, 'rt');
  
  if clean == -1 || raw == -1 || count == -1
      error(['Can''t open csv file for ' prefix ' in ' directory]);
  end
  
  linec = fgetl(clean);
  liner = fgetl(raw);
  linet = fgetl(count);
  
  linec = trim(linec);
  liner = trim(liner);
  linet = trim(linet);
  
  if isempty(linec)
      fprintf ('clean: %s\n', cleanfile);
      fprintf ('raw:   %s\n', rawfile);
      fprintf ('count: %s\n', countfile);
      
      error('CSV file header missing')
  end
  
  if (~strcmp(linec, liner) || ~strcmp(linec, linet))
      fprintf ('clean: %s\n', cleanfile);
      fprintf ('raw:   %s\n', rawfile);
      fprintf ('count: %s\n', countfile);
      
     error 'Header line mismatch';
  end

  if (~isempty(linelast) && ~strcmp(linelast,linec))
      fprintf ('clean: %s\n', cleanfile);
      fprintf ('raw:   %s\n', rawfile);
      fprintf ('count: %s\n', countfile);
      
     error(['Header line differs from previous section ' sections{f} ])
  end

  % get columns from header
  
  fields = strtrim(regexp(liner, ',', 'split'));

  for i = 1:length(fields)
    if (strcmp(fields{i},'Interval'))
      intrval = i;
    end
    if (strcmp(fields{i},'Layer'))
      layer = i;
    end
    if (strcmp(fields{i},'Layer_depth_min'))
      layer_min = i;
    end
    if (strcmp(fields{i},'Layer_depth_max'))
      layer_max = i;
    end
    if (strcmp(fields{i},'Samples'))
       samples = i;
    end
    if (strcmp(fields{i},'Good_samples'))
       samples = i;
    end
    if (strcmp(fields{i},'Lat_M'))
       latitude = i;
    end
    if (strcmp(fields{i},'Sv_mean'))
      sv_mean = i;
    end
    if (strcmp(fields{i},'Standard_deviation'))
      sv_sd = i;
    end
    if (strcmp(fields{i},'Skewness'))
      sv_skew = i;
    end
    if (strcmp(fields{i},'Kurtosis'))
      sv_kurt = i;
    end
  end
  
  if (samples == 0)
     error 'Samples column not found'
  end

  if (layer == 0)
     error 'Layer column not found'
  end
  
  if (first)
    extras = [ ', Linear_sv, Pct_good, Layer_depth, Uncleaned_Linear_sv' ...
        ', Uncleaned_Sv_mean, Uncleaned_Standard_deviation' ...
        ', Uncleaned_Skewness, Uncleaned_Kurtosis' ];
    fprintf(out, '%s%s\n',  linec, extras );
    first = 0;
  end

  cdata = textscan(clean,'%s','Delimiter','');
  cdata = nsort(cdata{1});
  rdata = textscan(raw,'%s','Delimiter','');
  rdata = nsort(rdata{1});
  tdata = textscan(count,'%s','Delimiter','');
  tdata = nsort(tdata{1});
  
  fclose(clean);
  fclose(raw);
  fclose(count);

  cleanline=0;
  rawline=0;
  countline=0;
  
  % process data
  while cleanline < length(cdata)
    cleanline = cleanline+1;
    rawline = rawline+1;
    countline = countline +1;
     
    linec = cdata{cleanline};
    liner = rdata{rawline};
    linet = tdata{countline};
    cfields = regexp(linec, ',', 'split');
    rfields = regexp(liner, ',', 'split');
    tfields = regexp(linet, ',', 'split');
    
    % skip lines with no position
    lat = strtrim(cfields{latitude});
    if strcmp(lat(1:3), '999')
      rawline = rawline-1;
      countline = countline-1;
      continue;
    end
      
    cinterval = str2double(cfields{intrval});
    clayer = str2double(cfields{layer});
    
    rinterval = str2double(rfields{intrval});
    rlayer = str2double(rfields{layer});
        
    % skip raw data without clean data
    while rawline < length(rdata) && ...
        (cinterval > rinterval || clayer > rlayer)
        rawline = rawline+1;
        rfields = regexp(rdata{rawline}, ',', 'split');
        rinterval = str2double(rfields{intrval});
        rlayer = str2double(rfields{layer});
    end        
    
    tinterval = str2double(tfields{intrval});
    tlayer = str2double(tfields{layer});
    
    % skip count data without clean data
    while countline < length(tdata) &&  ...
        (cinterval > tinterval || clayer > tlayer)
        countline = countline +1;
        tfields = regexp(tdata{countline}, ',', 'split');
        tinterval = str2double(tfields{intrval});
        tlayer = str2double(tfields{layer});
    end
    
    if cinterval ~= rinterval || clayer ~= rlayer || ...
            cinterval ~= tinterval || clayer ~= tlayer
        fprintf ('clean: %s\n', cleanfile);
        fprintf ('raw:   %s\n', rawfile);
        fprintf ('count: %s\n', countfile);
 
        error('file synchronisation lost at interval: %d %d %d layer: %d %d %d', ...
            cinterval, rinterval, tinterval, clayer, rlayer, tlayer);
    end
    
    % percent good
    cleandata = str2double(tfields{samples});
    if (cleandata > 0)
        rawdata = str2double(rfields{samples});
        pctgood = floor(100 * cleandata / rawdata);
    else
        pctgood = 0;
    end
    
    % layer depth
    if (layer_min > 0 && layer_max > 0)
        layer_depth = (str2double(cfields{layer_min}) + ...
            str2double(cfields{layer_max})) / 2;
    else
        layer_depth = clayer * 10 - 5;
    end
    
    % skip data which doesn't satisfy threshold conditions
    if (pctgood < MIN_GOOD) || (layer_depth > MAX_DEPTH)
        continue;
    end
    
    % convert to linear
    csv = str2double(cfields{sv_mean});
    if csv == 0
        csv = 9999;
    end
    if csv ~= 9999
        csv = 10 ^ (csv / 10);
    end
    
    rsv = str2double(rfields{sv_mean});
    if rsv == 0
        rsv = 9999;
    end
    if rsv ~= 9999
        rsv = 10 ^ (rsv / 10);
    end
    
    % generate output line
    lineout = [ linec ', ' num2str(csv) ', ' int2str(pctgood) ', ' ...
        num2str(layer_depth) ', ' num2str(rsv)...
        ',' rfields{sv_mean} ',' rfields{sv_sd} ',' rfields{sv_skew} ',' rfields{sv_kurt}];

    % if a new interval output the last interval
    if (cinterval > pending)
        for l = 1:currlayer
            if ~isempty(buffer{l})
                fprintf(out, '%s\n', buffer{l});
                buffer{l} = [];
            end
        end
      pending = cinterval;
      currlayer = 0;
    end

    % if current interval buffer (or overwrite buffer) this layer
    if (cinterval == pending)
        currlayer = clayer;
        buffer{currlayer} = lineout;
    end

  end

end

% output last buffer
for l = 1:currlayer
    if ~isempty(buffer{l})
        fprintf(out, '%s\n', buffer{l});
        buffer{l} = [];
    end
end

fclose(out);

end

function line = trim(line)
% Remove UTF-8 Byte Order Mark (ef bb bf) from header line if present
%

if  strncmp(line, ['' 239 187 191], 3)          % UTF-8 - used by EchoView 4 & 5
    line = line(4:end);
end

line=strtrim(line);

end

function lines = nsort(lines)
%sort comma separated lines by numeric values in fields

index = zeros(length(lines),3);
for i = 1 : length(lines)
    index(i,3) = i;
    line = lines{i};
    c = find(line == ',',3);
    index(i,1) = str2double(line(c(1)+1:c(2)-1));
    index(i,2) = str2double(line(c(2)+1:c(3)-1));
end
index = sortrows(index);
lines = lines(index(:,3));
end
