% GET_SYNTS_PROFILES  Quick and dirty (but not fast) way to extract synTS
%    profiles at a set of locations and depths, for one target date.
%
% INPUT:
%   tim     - single target time (decimal days since 1900), or [year month day]
%   lon,lat - vector of N locations
%   deps    - vector of D depths
%   tdel    - will search for synTS file at tim+/-tdel (days) [default 2]
%   conf    - structure with optional fields:
%     get_t - 1=get temperature   [default 1]
%     get_s - 1=get salinity      [default 1]
%
% OUTPUT
%   tt      - [D N] temperature profiles  (empty if get_t=0)
%   ss      - [D N] salinity profiles  (empty if get_t=0)
%   age     - synTS estimate is centred on 'age' days from tim
%
% AUTHOR:  Jeff Dunn  CMAR   28/2/08
%
% USAGE: [tt,ss,age] = get_synTS_profiles(tim,lon,lat,deps,tdel,conf);

function [tt,ss,age] = get_synTS_profiles(tim,lon,lat,deps,tdel,conf)

tt = [];
ss = [];

if length(tim)==3
   tim = greg2time([tim 0 0 0]);
end

if nargin<5 || isempty(tdel)
   tdel = 7;
end
if nargin<6 || isempty(conf) || ~isfield(conf,'get_t') || isempty(conf.get_t)
   conf.get_t = 1;
end
if ~isfield(conf,'get_s') || isempty(conf.get_s)
   conf.get_s = 1;
end

   

pth{1} = '/home/datalib/climatologies/synTS/hindcast13/';
pth{2} = '/home/datalib/climatologies/synTS/hindcast06/';
pth{3} = '/home/datalib/climatologies/synTS/NRT06/';
pth{4} = '/home/datalib/climatologies/synTS/NRT06_v1/';

toff(2*(1:tdel)) = -(1:tdel);
toff(2*(1:tdel)+1) = (1:tdel);


ipth = 1;
ii = 1;
found = 0;

while ipth<=length(pth) && ~found
   dstr = time2greg(tim+toff(ii));
   fnm = sprintf('%ssynTS_%4d%0.2d%0.2d.nc',pth{ipth},dstr(1:3));
   found = exist(fnm,'file'); 
		  
   if ~found
      ii = ii+1;
      if ii>length(toff)
	 ii = 1;
	 ipth = ipth+1;
      end
   end
end		  
age = toff(ii);


if found
   x = getnc(fnm,'lon');
   y = getnc(fnm,'lat');
   dep = getnc(fnm,'depth');
   
   x1 = max(find(x<min(lon)));
   if isempty(x1); x1 = 1; end
   x2 = min(find(x>max(lon)));
   if isempty(x2); x2 = length(x); end
   y1 = max(find(y<min(lat)));
   if isempty(y1); y1 = 1; end
   y2 = min(find(y>max(lat)));
   if isempty(y2); y2 = length(y); end
   
   [x,y,z] = meshgrid(x(x1:x2),y(y1:y2),dep);
   [X,Z] = meshgrid(lon,deps);
   [Y,Z] = meshgrid(lat,deps);
   
   if conf.get_s
      s = shiftdim(squeeze(getnc(fnm,'salinity',[-1 -1 y1 x1],[-1 -1 y2 x2])),1);
      ss = interp3(x,y,z,s,X,Y,Z);
      clear s
   end
   
   if conf.get_t
      t = shiftdim(squeeze(getnc(fnm,'temperature',[-1 -1 y1 x1],[-1 -1 y2 x2])),1);
      tt = interp3(x,y,z,t,X,Y,Z);
      clear t
   end
end

%------------------------------------------------------------------------------
