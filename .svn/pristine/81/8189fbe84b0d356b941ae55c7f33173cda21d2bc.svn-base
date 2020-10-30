function sample_data = get_woa98(sample_data, depth_offset, tfile, sfile)

% WRITTEN BUT NOT YET WORKING - some debugging still to do.

% ftp://ftp.cdc.noaa.gov/Datasets/nodc.woa98/temperat/monthly/otemp.anal1deg.nc
% ftp://ftp.cdc.noaa.gov/Datasets/nodc.woa98/salinity/monthly/salt.anal1deg.nc

% get_woa98 Get temperature and salinity data from the World Ocean Atlas 98
% (or compatible) dataset and add it an IMOS-toolbox sample_data structure.
% 
% sample_data must have:
% a TIME and a DEPTH dimension,
% LATITUDE and LONGITUDE variables of TIME dimension,
%
% depth_offset is added to the sample_data DEPTH to get the ocean depth
% (sample_data is relative to transducer, woa98 is relative to surface).
% 
% The output sample_data will have temperature and salinity variables added.
% temperature and salinity haved dimensions TIME, DEPTH.
% These are profiles interpolated from synTS data in space and time to the 
% position given by TIME, LATITUDE, LONGITUDE and DEPTH (+ depth_offset).

if nargin < 2
    depth_offset = 0;
end

% Get Temperature and salinity files
if nargin < 3
    tfile = 'otemp.anal1deg.nc';
end
if nargin < 4
    sfile = 'salt.anal1deg.nc';
end
if exist(tfile,'file') ~= 2
    [file,path] = uigetfile('*.nc','WOA 98 Temperature data',tfile);
    tfile = fullfile(path,file);
end
if exist(sfile, 'file') ~= 2
    [file,path] = uigetfile('*.nc','WOA 98 Salinity data', ...
        fullfile(fileparts(tfile),sfile));
    sfile = fullfile(path,file);
end

% identify time, latitude and longitude in sample_data
timed = 0;
depthd = 0;
latv = 0;
lonv = 0;
tempv = 0;
salv = 0;
for k = 1:length(sample_data.dimensions)
    if strcmp(sample_data.dimensions{k}.name, 'TIME')
        timed = k;
    end
    if strcmp(sample_data.dimensions{k}.name, 'DEPTH')
        depthd = k;
    end
end
for k = 1:length(sample_data.variables)
    if strcmp(sample_data.variables{k}.name, 'LATITUDE')
        latv = k;
    end
    if strcmp(sample_data.variables{k}.name, 'LONGITUDE')
        lonv = k;
    end
    if strcmp(sample_data.variables{k}.name, 'temperature')
        tempv = k;
    end
    if strcmp(sample_data.variables{k}.name, 'salinity')
        salv = k;
    end
end

stime = sample_data.dimensions{timed}.data;
sdepth = double(sample_data.dimensions{depthd}.data) + depth_offset;
slat = sample_data.variables{latv}.data;
slon = sample_data.variables{lonv}.data;
slon(slon < 0) = slon(slon < 0) + 360;

mnd = min(sdepth);
mxd = max(sdepth);
minlat = min(slat);
maxlat = max(slat);
minlon = min(slon);    
maxlon = max(slon);

% Note: if data crosses longitude 0 we will extract more data than we need
% (all longitudes) but the code should still work. Added support for 0
% crossings is therefore an efficiency enhancement not considered
% worthwhile at this point in time.


% Read netcdf data
tncid = netcdf.open(tfile, 'NC_NOWRITE');

tid = netcdf.inqVarID(tncid, 'time');
did = netcdf.inqVarID(tncid, 'level');
latid = netcdf.inqVarID(tncid, 'lat');
lonid = netcdf.inqVarID(tncid, 'lon');

nctime = double(netcdf.getVar(tncid, tid));
ncdepth = double(netcdf.getVar(tncid, did));
nclat = double(netcdf.getVar(tncid, latid));
nclon = double(netcdf.getVar(tncid, lonid));

% assumption - temp and salinity files have the same lat, lon, level, time

% start, finish and count are 0 indexed for netcdf
[year, month, ~] = datevec(stime);
sday = stime - datenum(year,1,1);

start = [ max(1,find(nclon > minlon,1,'first')-2) ...
    max(1,find(nclat < maxlat,1,'first')-2) ...
    max(1,find(ncdepth > mnd,1,'first')-2) ...
    min(month) - 1];
finish = [ min(length(nclon), find(nclon < maxlon,1,'last')) ...
    min(length(nclat), find(nclat > minlat,1,'last')) ...
    min(length(ncdepth), find(ncdepth < mxd,1,'last')) ...
    max(month)];
count = finish - start + 1;

lonx = 1 + (start(1):finish(1));
latx = 1 + (start(2):finish(2));
dpx = 1 + (start(3):finish(3));
tmx= 1 + (start(4):finish(4));

% read data of interest from netcdf file
sncid = netcdf.open(sfile, 'NC_NOWRITE');
tempid = netcdf.inqVarID(tncid, 'otemp');
saltid = netcdf.inqVarID(sncid, 'salt');

if tmx(end) == 13
    cnt = count;
    cnt(4) = cnt(4) - 1;
    strt = start;
    strt(4) = 0;
    cnt1 = count;
    cnt1(4) = 1;
    nctemp = [double(netcdf.getVar(tncid, tempid, start, cnt)) ...
        double(netcdf.getVar(tncid, tempid, strt, cnt1))];
    ncsalt = [double(netcdf.getVar(sncid, saltid, start, cnt)) ...
        double(netcdf.getVar(sncid, saltid, strt, cnt1))];
else
    nctemp = double(netcdf.getVar(tncid, tempid, start, count));
    ncsalt = double(netcdf.getVar(sncid, saltid, start, count));
end
netcdf.close(tncid);
netcdf.close(sncid);

% Interpolate temp and salt at each depth level for each location.

itemp = nan(length(stime),length(dpx));
isalt = nan(length(stime),length(dpx));
for i = dpx
    itemp(:,i) = interp3(nclon(lonx),nclat(latx),nctime(tmx), squeeze(nctemp(:,:,i,:)),slon,slat,sday); %TODO
    isalt(:,i) = interp3(nclon(lonx),nclat(latx),nctime(tmx), squeeze(ncsalt(:,:,i,:)),slon,slat,sday); %TODO
end

% Interpolate temp and salt at sdepth depths

temperature = interp2(stime,ncdepth,itemp,stime,sdepth);
salinity = interp2(stime,ncdepth,isalt,stime,sdepth);

% store results in sample_data. If variable does not exist create it.
if tempv == 0
    tempv = length(sample_data.variables)+1;
    sample_data.variables{tempv}.name = 'temperature';
    sample_data.variables{tempv}.dimensions(1) = timed;
    sample_data.variables{tempv}.dimensions(2) = depthd;
end
sample_data.variables{tempv}.source = 'WOA98';
sample_data.variables{tempv}.data = temperature';

if salv == 0
    salv = length(sample_data.variables)+1;
    sample_data.variables{salv}.name = 'salinity';
    sample_data.variables{salv}.dimensions(1) = timed;
    sample_data.variables{salv}.dimensions(2) = depthd;
end
sample_data.variables{salv}.source = 'WOA98';
sample_data.variables{salv}.data = salinity';


% update history
nowj = (now - datenum([1970 1 1])) * 86400000;              % now in ms since 1970
timezone = java.util.TimeZone.getDefault().getOffset(nowj); % timezone offset in ms  
nowt = now - timezone / 86400000;                           % now UTC in days       
comment = [datestr(nowt, 'yyyy-mm-ddTHH:MM:SSZ') ' '...
    getenv('USER') getenv('UserName') ...
    ' Inferred temperature and salinity read from WOA 98'];
if isfield(sample_data, 'history') && ~isempty(sample_data.history);
    comment = sprintf('%s\n%s', sample_data.history, comment);
end
sample_data.date_modified = nowt;
sample_data.history = comment;




