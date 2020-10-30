function sample_data = get_clim(sample_data, depth_offset, woafile)
% Get temperate and salinity from woa98Climatology.nc, Lindsay Pender's
% combination of World Ocean Atlas 98 and CARS 2000.
%
% Inputs:
%       sample_data     IMOS toolbox format sample_data structure which
%                       must have TIME and DEPTH dimensions and LATITUDE
%                       and LONGITUDE variables of dimension TIME.
%       depth_offset    Amount to add to sample_data depth to get real
%                       depth (transducer depth) [0]
%       woafile         Full path to woa98Climatology.nc. If not found will
%                       ask the user unless a noncharacter value is given.
%                       [woa98Climatology.nc]
%
% Outputs:
%       sample_data     sample_data with temperature and salinity variables
%                       with data from woafile. sample_data.history is also
%                       updated.
%
% Author: Gordon Keith
% Date: 20141029

ncfile = 'woa98Climatology.nc';

if nargin < 2
    depth_offset = 0;
end

if nargin < 3
    woafile = ncfile;
else
    if ischar(woafile) && ~isempty(woafile)
        ncfile = woafile;
    end
end

if ~exist(ncfile, 'file')
    if ~ischar(woafile)
        error('World Oceanography Atlas not found')
    end
    [file,path] = uigetfile('*.nc','Locate WOA file',ncfile);
    if isequal(file, 0)
        error('World Oceanography Atlas not found')
    end
    ncfile = fullfile(path,file);
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


% Open ncfile
ncid = netcdf.open(ncfile, 'NC_NOWRITE');

did = netcdf.inqVarID(ncid, 'depth');
latid = netcdf.inqVarID(ncid, 'latitude');
lonid = netcdf.inqVarID(ncid, 'longitude');


ncdepth = double(netcdf.getVar(ncid, did));
nclat = double(netcdf.getVar(ncid, latid));
nclon = double(netcdf.getVar(ncid, lonid));

% assumption - latitude and longitude are regular arrays
lat0 = nclat(1);
dlat = nclat(2) - nclat(1);
lon0 = nclon(1);
dlon = nclon(2) - nclon(1);

% start, finish and count are 0 indexed for netcdf
start = [ max(1,find(nclon > minlon,1,'first')-2) ...
    max(1,find(nclat > minlat,1,'first')-2) ...
    max(1,find(ncdepth > mnd,1,'first')-2) ];
finish = [ min(length(nclon), find(nclon < maxlon,1,'last')) ...
    min(length(nclat), find(nclat < maxlat,1,'last')) ...
    min(length(ncdepth), find(ncdepth < mxd,1,'last'))];
dpx = start(3):finish(3);
dx = 1:length(dpx);
count = finish - start + 1;

% read data of interest from netcdf file
mt = getVar(ncid, 'temperature', start, count);
ast = getVar(ncid, 'anSinTemperature', start, count);
act = getVar(ncid, 'anCosTemperature', start, count);
sst = getVar(ncid, 'saSinTemperature', start, count);
sct = getVar(ncid, 'saCosTemperature', start, count);
ms = getVar(ncid, 'salinity', start, count);
ass = getVar(ncid, 'anSinSalinity', start, count);
acs = getVar(ncid, 'anCosSalinity', start, count);
sss = getVar(ncid, 'saSinSalinity', start, count);
scs = getVar(ncid, 'saCosSalinity', start, count);
netcdf.close(ncid);

% calculate temperature and salinity for each profile
temperature = nan(length(stime),length(sdepth));
salinity = nan(length(stime),length(sdepth));
do=ones(size(sdepth));

[year,~,~] = datevec(stime(1));
st = (stime - datenum(year,1,1)) / 365 * 2 * pi; % date in radians of year

for i = 1:length(stime)
    % select box surrounding location
    latx = floor((slat(i) - lat0) / dlat);
    lonx = floor((slon(i) - lon0) / dlon);
    latx = [latx latx+1];                   %#ok<AGROW>
    lonx = [lonx lonx+1];                   %#ok<AGROW>
    lnx = lonx - start(1) + 1;
    ltx = latx - start(2) + 1;
    
    % calculate temperature and salinity for surrounding box then
    % interpolate at profile location.
    tmp = mt(lnx,ltx,dx) + sin(st(i)) * ast(lnx,ltx,dx) + cos(st(i)) * act(lnx,ltx,dx) + ...
        sin(2*st(i)) * sst(lnx,ltx,dx) + cos(2*st(i)) * sct(lnx,ltx,dx);
    temperature(i,:) = interp3(nclon(lonx),nclat(latx),ncdepth(dpx), ...
        tmp, slon(i)*do, slat(i)*do, sdepth);
    
    slt = ms(lnx,ltx,dx) + sin(st(i)) * ass(lnx,ltx,dx) + cos(st(i)) * acs(lnx,ltx,dx) + ...
        sin(2*st(i)) * sss(lnx,ltx,dx) + cos(2*st(i)) * scs(lnx,ltx,dx);
    salinity(i,:) = interp3(nclon(lonx),nclat(latx),ncdepth(dpx), ...
        slt, slon(i)*do, slat(i)*do, sdepth);     
end

% store results in sample_data. If variable does not exist create it.
if tempv == 0
    tempv = length(sample_data.variables)+1;
    sample_data.variables{tempv}.name = 'temperature';
    sample_data.variables{tempv}.dimensions(1) = timed;
    sample_data.variables{tempv}.dimensions(2) = depthd;
end
sample_data.variables{tempv}.source = ncfile;
sample_data.variables{tempv}.data = temperature';

if salv == 0
    salv = length(sample_data.variables)+1;
    sample_data.variables{salv}.name = 'salinity';
    sample_data.variables{salv}.dimensions(1) = timed;
    sample_data.variables{salv}.dimensions(2) = depthd;
end
sample_data.variables{salv}.source = ncfile;
sample_data.variables{salv}.data = salinity';

% update history
nowj = (now - datenum([1970 1 1])) * 86400000;              % now in ms since 1970
timezone = java.util.TimeZone.getDefault().getOffset(nowj); % timezone offset in ms  
nowt = now - timezone / 86400000;                           % now UTC in days       
comment = [datestr(nowt, 'yyyy-mm-ddTHH:MM:SSZ') ' '...
    getenv('USER') getenv('UserName') ...
    ' Inferred temperature and salinity read from World Ocean Atlas 98 and CARS 2000'];
if isfield(sample_data, 'history') && ~isempty(sample_data.history);
    comment = sprintf('%s\n%s', sample_data.history, comment);
end
sample_data.date_modified = nowt;
sample_data.history = comment;

function var = getVar(ncid, name, start, count)
v_id = netcdf.inqVarID(ncid, name);
var = double(netcdf.getVar(ncid, v_id, start, count));
fill = netcdf.getAtt(ncid,v_id,'_FillValue');
offset = netcdf.getAtt(ncid,v_id,'add_offset');
scale = netcdf.getAtt(ncid,v_id,'scale_factor');
var(var == fill) = nan;
var = offset + scale * var;



