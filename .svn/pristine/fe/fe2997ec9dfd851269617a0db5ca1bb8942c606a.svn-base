% quick_compo_Wsat:  Construct a quick-and-dirty composite of Windsat SST v7,
% by just filling gaps from most recent daily file with data from next most
% recent, and so on for ndays
%
% INPUT: tdate   Target date, days since 1900
%        ndays   Total number of days data to use
%        lo,la   Locations of required SST
%
% OUTPUT composite SST
%
% Jeff Dunn 26/3/12
%
% USAGE: sso = quick_compo_Wsat(tdate,ndays,lo,la);

function sso = quick_compo_Wsat(tdate,ndays,lo,la)

% MODS:  11/6/2013 JRD Set path according to year (hardcoded to 2012 before)

sso = [];

fpth = '/home/datalib/platforms/wsat/unzipped/';


for ii = 1:ndays
   dat = time2greg(tdate-(ii-1));
   syr = sprintf('%4d',dat(1));
   sdate = sprintf('%4d%0.2d%0.2d',dat(1:3));
   fnm = [fpth syr '/' sdate '-wsat-remss-l2p_gridded_25-wsat_' sdate 'v7-v01.nc'];
   if ~exist(fnm,'file')
      fnm = [fpth syr '/' sdate '-wsat-remss-l2p_gridded_25-wsat_' sdate 'rt-v01.nc'];
      if ~exist(fnm,'file')
	 fnm = [];
      end
   end
   
   if ~isempty(fnm)
      ncload(fnm,'sea_surface_temperature','confidence_flag','proximity_confidence');
      % USe first (descending) passes as these are early morning, hence not
      % affected by skin heating 
      ss = squeeze(sea_surface_temperature(1,:,:));
      ss(ss<-32767)=nan;
      ss(squeeze(confidence_flag(1,:,:))~=0) = nan;
      ss(squeeze(proximity_confidence(1,:,:))<=2) = nan;
      
      if isempty(sso)
	 sso = ss./100;
	 ncload(fnm,'lon','lat');
	 [x,y] = meshgrid(lon,lat);
      else
	 jj = isnan(sso);
	 sso(jj) = ss(jj)./100;
      end
   end
end
   
if isempty(sso)
   sso = nan(size(la));
else
   sso = spread1cell(sso);      % dirty - spread data edges to improve interp recovery
   sso = interp2(x,y,sso,lo,la);
end
