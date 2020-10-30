function [sst,lon_used,lat_used,dates] = pftseries(pos,varopt)
% Retrieves a time_series of 10-day averaged estimates 
% of NOAA/NASA Pathfinder SST
% closest to the lat/lon requested.
%
%  [sst,lon_used,lat_used,dates] = ptseries([lon lat],varopt);
%-----------------------------------------------------------------------
% POS [lon lat] :is a 2-element vector of longitude and latitude
%
% VAROPT can be:
%         'sst' to get the SST estimate (default)
%         'err' to get the expected error
%         'raw' to get the unprocessed Pathfinder best-SST
%
% LON_used, LAT_used: actual location of estimates
%
% DATES: days since 01/01/1985
%
% SST:	 sea surface temperature in degrees celcius
%----------------------------------------------------------------------
% This version is configured for a dataset with boundaries 0-70 S and
% 90-200 E and posible date ranges are between 24 Feb 1987 and June 1994
%
% Alison Walker Aug 97
%---------------------------------------------------------------------------

% CHECK INPUTS

if nargin ==2 
  if ~(strcmp(varopt,'sst') | strcmp(varopt,'err') | strcmp(varopt,'raw'))
    error(['varopt must be sst, err or raw.  It was: ' varopt])
  end
else
  varopt = 'sst'; % default is to get the sst estimate
end

% this logical value gets tested repeatedly
reading_raw_data = strcmp(varopt,'raw');

lon = pos(1);
lat = pos(2);

if lon < 90
   error(['Check limits: min longitude < ' num2str(90)])
elseif lon > 200
  error(['Check limits: max longitude > 200'])

elseif lat > 0
  error(['Check limits: max latitude > 0'])
elseif lat < -70
  error(['Check limits: min latitude < -70'])
end

% FIND LAT, LON INDICES

% 1/delta_angle of data grid interval
delta=4096./360.;
% index limits
ypix = -round(lat*delta)+1;
xpix = round((lon-90)*delta)+1;

% CHOOSE APPROPRIATE NETCDF FILE

year = 87;
if reading_raw_data
   prefix=(['/home/eez_data/sst/pathfinder_aus/pfsstr']);
else
   prefix=(['/home/eez_data/sst/pathfinder_aus/pfsste']);
end
fname=([prefix num2str(year)]);



% RETRIEVE REQUESTED DATA FROM NETCDF FILES

% retrieve lat/lon data from netcdf file
lat_used=getnc(fname,'lat',[ypix],[ypix],[1]);
lon_used=getnc(fname,'lon',[xpix],[xpix],[1]);

% retrieve sst or error data from netcdf files
sst = []; dates = [];
for yr = 87:94
   fname=([prefix num2str(yr)]);
   sub_sst = [];
   sub_dates = getnc(fname,'time');
   len = length(sub_dates);
   if strcmp(varopt,'err')
  	sub_sst =getnc(fname,'error',...
	 [1 ypix xpix],[len ypix xpix],[1 1 1]);
   else
   	sub_sst =getnc(fname,'sst',...
	 [1 ypix xpix],[len ypix xpix],[1 1 1]);
   end
   sst = [sst; sub_sst];
   dates = [dates; sub_dates];
end

if reading_raw_data
   sst=.15*sst-3.0;
   sst=change(sst,'==',-3.0,nan);
end
