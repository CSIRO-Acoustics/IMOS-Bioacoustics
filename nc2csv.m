function data = nc2csv(ncfile, csvfile, channel)
% Read an IMOS BASOOP netcdf file and write an echoview .sv.csv file of the
% Sv data.
%
% Inputs
%   ncfile  full path to netcdf file [will ask if not provided]
%   csvfile name of csv file to create [will ask if not provided]
%   channel Which channel to write [1]
%
% Author Gordon Keith 20151020

if nargin < 1
    ncfile = '';
end

if exist(ncfile, 'file') ~= 2
    [file, path] = uigetfile('*.nc', 'IMOS BASOOP netcdf file', ncfile);
    if isequal(file,0)
        error('NetCDF file required');
    end
    ncfile = fullfile(path,file);
end

if nargin < 2 || isempty(csvfile)
    csvfile = fileparts(ncfile);
end
if exist(csvfile, 'dir') == 7
    [file, path] = uiputfile('*.sv.csv', 'Sv csv file to save', fileparts(ncfile));
    if isequal(file,0)
        error('.sv.csv file name required');
    end
    if ~strcmpi(file(end-3:end), '.csv')
        file = [file '.sv.csv'];
    end
    csvfile = fullfile(path,file);
end

if nargin < 3
    channel = 1;
end

data = viz_sv(ncfile,[],'noplot');

gfile = csvfile;
if strcmpi(gfile(end-3:end), '.csv')
    gfile(end-3:end) = [];
end
if strcmpi(gfile(end-2:end), '.sv')
    gfile(end-2:end) = [];
end
gfile = [gfile '.gps.csv'];

fid = fopen(csvfile,'w');
gid = fopen(gfile, 'w');

if fid < 0
    error('Unable to open %s', csvfile);
end

fprintf(fid, 'Ping_date, Ping_time, Range_start, Range_stop, Sample_count,\n');
fprintf(gid, 'GPS_date, GPS_time, Latitude, Longitude\n');

ns = length(data.depth);
h = (data.depth(end) - data.depth(1)) / (ns-1);
range = sprintf('%g, %g, %g', data.depth(1) - h/2, data.depth(end) + h/2, ns);

for i = 1:length(data.time)
    dt = datestr(data.time(i), 'yyyy-mm-dd, HH:MM:SS');
    fprintf(fid, '%s, %s', dt, range);
    fprintf(fid, ', %g', data.Sv(:,i,channel));
    fprintf(fid, '\n');
    fprintf(gid, '%s, %g, %g\n', dt, data.latitude(i), data.longitude(i));
end

fclose(fid);
fclose(gid);
