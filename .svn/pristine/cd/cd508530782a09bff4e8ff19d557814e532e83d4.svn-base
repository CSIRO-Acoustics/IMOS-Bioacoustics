function sample_data = get_npp(sample_data,npp_path,days)
% Add Net Primary Production variable to IMOS sample_data structure
%
% Inputs:
%   sample_data - IMOs sample data structure 
%   npp_path - directory holding year/npp.*.hdf files
%   days - number of days to integrate npp over [365]



if nargin < 2
    settings = basoop();
    npp_path = settings.npp_path;

end

if nargin < 3
    days = 365;
end


[year, ~] = datevec(sample_data.time_coverage_end);

endtime = (sample_data.time_coverage_end - datenum(1970,1,1)) * 86400;
yearstart = endtime - days * 86400;

fy = dir(fullfile(npp_path,num2str(year),'npp.*.hdf'));
fy_1 = dir(fullfile(npp_path,num2str(year - 1),'npp.*.hdf'));

% Edit by Haris to read 'vgpm.*.hdf' files; 10 May 2018

if isempty(fy)
    fy = dir(fullfile(npp_path,num2str(year),'vgpm.*.hdf')); 
end

if isempty(fy_1)
    fy_1 = dir(fullfile(npp_path,num2str(year - 1),'vgpm.*.hdf'));
end

files = {};
x = -1;
y = -1;
info = [];

for file = 1:length(fy_1)
    ffile = fullfile(npp_path,num2str(year - 1),fy_1(file).name);
    inf=hdfinfo(ffile);
    if length(inf.Attributes) > 1 && ...
        strcmp(inf.Attributes(1).Name, 'Start Time') && ...
        inf.Attributes(1).Value < endtime && ...
        inf.Attributes(1).Value > yearstart
        if x < 0 
            x = inf.SDS.Dims(2).Size;
            y = inf.SDS.Dims(1).Size;
            info = inf;
        elseif  x ~= inf.SDS.Dims(2).Size || ...
                y ~= inf.SDS.Dims(1).Size;
            error('File sizes do not match');
        end
        files{end+1}=ffile;
    end
end

for file = 1:length(fy)
    ffile = fullfile(npp_path,num2str(year),fy(file).name);
    inf=hdfinfo(ffile);
    if length(inf.Attributes) > 1 && ...
        strcmp(inf.Attributes(1).Name, 'Start Time') && ...
        inf.Attributes(1).Value < endtime && ...
        inf.Attributes(1).Value > yearstart
        if x < 0 
            x = inf.SDS(1).Dims(2).Size;
            y = inf.SDS(1).Dims(1).Size;
            info = inf;
        elseif  x ~= inf.SDS(1).Dims(2).Size || ...
                y ~= inf.SDS(1).Dims(1).Size;
            error('File sizes do not match');
        end
        files{end+1}=ffile;
    end
end

if isempty(files)
    error('Could not find any suitable files');
end

limit = info.SDS(1).Attributes(7).Value;

if ~strcmp(info.SDS(1).Attributes(7).Name, 'Limit') 
    warning('File assumptions don''t match');
    limit = [-90 -180 90 180];
end

% NPP data has image orientation, row 1 is north, row end is south 

west = floor((sample_data.geospatial_lon_min - limit(2)) * x / (limit(4) - limit(2)));
east = ceil ((sample_data.geospatial_lon_max - limit(2)) * x / (limit(4) - limit(2)));
south = ceil((-sample_data.geospatial_lat_min - limit(1)) * y / (limit(3) - limit(1)));
north = floor((-sample_data.geospatial_lat_max - limit(1)) * y / (limit(3) - limit(1)));

bound_adj = 0;
if east < west
    bound_adj = x;
end

npp = nan(length(files), south - north + 1, east - west + bound_adj + 1);
cnt = zeros(size(npp));

for file = 1:length(files)
    inf = hdfinfo(files{file});
    days = (inf.Attributes(2).Value - inf.Attributes(1).Value ) / 86400 + 1;
    if bound_adj == 0
        npp(file,:,:) = hdfread(files{file}, info.SDS(1).Name, 'Index', {[north west], [1 1], [south - north + 1 east - west + 1]});
    else    % longitude wraps around file boundary (e.g. crosses dateline)
        npp(file,:,1:bound_adj - west) = hdfread(files{file}, info.SDS(1).Name, 'Index', {[north west], [1 1], [south - north + 1 bound_adj - west ]});
        npp(file,:,bound_adj - west + 1 : end) = hdfread(files{file}, info.SDS(1).Name, 'Index', {[north 1], [1 1], [south - north + 1 east + 1]});
    end
    cnt(file,:,:) = ones(size(cnt(file,:,:))) * double(days);
end

cnt(npp == info.SDS(1).Attributes(16).Value) = 0;
npp(npp == info.SDS(1).Attributes(16).Value) = 0;
anpp = squeeze(sum(npp .* cnt,1) ./ sum(cnt,1));

dlon = (west:east + bound_adj)   * (limit(4) - limit(2)) / x + limit(2);
dlat = -(north:south) * (limit(3) - limit(1)) / y - limit(1); 
lon = [];
lat = [];

for i = 1:length(sample_data.variables)
    if strcmpi(sample_data.variables{i}.name, 'LATITUDE') % usually i=1
        latv = i;
        lat = sample_data.variables{i}.data;
        break;
    end
end

for i = 1:length(sample_data.variables)
    if strcmpi(sample_data.variables{i}.name, 'LONGITUDE') % usually i=2
        lon = sample_data.variables{i}.data;
        if bound_adj
            adj = lon < west * (limit(4) - limit(2)) / bound_adj + limit(2) - 2;
            lon(adj) = lon(adj) + 360;
        end
        break;
    end
end

nppi = interp2(dlon,dlat,anpp,lon,lat);

nppv = length(sample_data.variables) + 1;
for i = 1:length(sample_data.variables)
    if strcmpi(sample_data.variables{i}.name, 'npp')
        nppv = i;
        break;
    end
end

sample_data.variables{nppv}.name = 'npp';
sample_data.variables{nppv}.dimensions = sample_data.variables{latv}.dimensions;
sample_data.variables{nppv}.data = nppi;
sample_data.variables{nppv}.long_name = 'net_primary_production';
sample_data.variables{nppv}.units = 'mg C m-2 day-1';
sample_data.variables{nppv}.FillValue_ = -9999;

%% Haris - 18/03/2020. Independent checking of Gordon's code
%{

% Independent checking of Gordon's code using 'ltln2val' function available
% in Matlab mapping toolbox:
% https://au.mathworks.com/help/map/ref/ltln2val.html . This function can
% extract data grid values for specified locations, exactly similar to what
% is implemented above. This checking was required to verify Gordon's
% modification impelemented on 17/03/2020.

for file = 1:length(files)
    npp_check{file} = hdfread(files{file}, '/npp', 'Index', {[1  1],[1  1],[1080  2160]}); % low resolution file - usual case
end

allData = cat(3,npp_check{:});
allData(allData == info.SDS(1).Attributes(16).Value) = NaN; % fill value to NaN

npp_check_mean = nanmean(allData,3); % average NPP for the whole year

% extract data grid values for specified locations

R = georasterref('RasterSize', [1080  2160], ...
      'RasterInterpretation', 'cells', 'ColumnsStartFrom', 'north', ...
      'LatitudeLimits', [-90 90], 'LongitudeLimits', [-180 180]);

  
for i = 1:length(sample_data.variables)
    if strcmpi(sample_data.variables{i}.name, 'LATITUDE') 
        latitu = sample_data.variables{i}.data;
    end
end

for i = 1:length(sample_data.variables)
    if strcmpi(sample_data.variables{i}.name, 'LONGITUDE') 
        longi = sample_data.variables{i}.data;
    end
end

nppi_check = ltln2val(npp_check_mean, R, latitu, longi,'linear'); % extract data grid values for specified locations

% compare result - minor difference in result could be due to different method used for averaging 

figure;
subplot(2,1,1)
plot(nppi,'DisplayName','Gordon result') % Gordon's result
box on; grid on; hold;
plot(nppi_check,'DisplayName','Mapping toolbox') % mapping toolbox result
ylabel ('NPP')
title('Comparision between two methods')
legend

subplot(2,1,2) 
plot(nppi-nppi_check)
box on; grid on;
ylabel ('Difference between methods')
xlabel ('Interval')
%}
