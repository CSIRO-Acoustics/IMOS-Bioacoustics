% GET_CLIM: Return extracts from CARS, at a set of locations and depths, and
%      optionally day-of-years
%      SUPERCEDES 'getsection' and can be used instead of 'get_clim_casts'
%      or 'getchunk' (both routines are FASTER but less general).
%      Note that 'get_clim_casts' cannot return other CARS variables.
%
% INPUT
%  prop    property ('t','s','o','si','n', OR 'p')
%  lon,lat  vector or array of NN locations
%  deps    vector of DD depths (in m). Need not be contiguous.
%  OPTIONAL:
%  doy     matrix day-of-year corresponding to locations
%  var     1 = return CARS property (mean or seasonal depending on 'fns' and
%          whether or not 'doy' is non-empty.)   If not = 1, then return
%          other climatology variables, (see GETMAP 'vars')
%  fnms    name or cell array of names of climatology files to
%          access (in the order specified). 
%      OR  'best' to get CARS2000 best possible resolution
%      OR  'best07' to get CARS2006 best possible resolution
%          Default:  'cars2000' [will REMAIN for backwards compatibility!]
%  woa     1=WOA05  2=WOA98   used outside of CARS region  [def 0]
%  fll     1=values at all depths outside of land areas (only with selected 
%          files) [def 0]
%  fpth    to get map files from a non-standard disc location
%  fns     0=mean only, 1=+annual, 2=+semi-ann [default: semi-ann if available]
%  config    Struct with fields for any extra options to override
%  -nneg   1=replace negative values with zero (esp for nutrients) 
%            [default 0 for temperature, otherwise 1] 
%
% OUTPUT
%  vv      psuedo-casts extracted from CARS, dimensioned [size(lon) DD]
%
% NOTE:  if fll~=1 then the interpolation used here will still provide some
%     crude inshore values from adjacent data points. These may result in 
%     temperature inversions etc which are not actually present in CARS. If
%     this is a problem, use fll=1 or use GET_CLIM_CASTS.
%
% AUTHOR: Jeff Dunn  CSIRO DMR  May 1998
% $Id: get_clim.m,v 1.1 2001/06/20 04:08:17 dun216 Exp dun216 $
% MODS:  see below
% CALLS:  intp3jd (getmap (clname))  csl_dep
%
% USAGE: vv = get_clim(prop,lon,lat,deps,doy,var,fnms,woa,fll,fpth,fns,config);

function vv = get_clim(prop,lon,lat,deps,doy,var,fnms,woa,fll,fpth,fns,config)

% MODS:
%  1/5/2001  Fixed range so that use WOA outside [30 200 -70 10] instead of
%            [100 200 -70 0]. Also, bug in WOA stuff if 2D x&y, and olny 1 depth.
%  20/6/01   Mod interp3_clim to extract the minimum required amount of CARS, to
%            decrease time and memory requirements.
%  18/7/07   Option to cap vv at 0 (replace -ve values)
%  30/5/08   Give access to WOA05 outside CARS region
%  6/8/10    Change default CARS version section

persistent GET_CLIM_def_warn

ncquiet;

if ~isa(lon,'double')
   lon = double(lon);
   lat = double(lat);
end

if nargin<5; doy = []; end
if length(doy)==1 & prod(size(lon))>1
   doy = repmat(doy,size(lon));
end
if nargin<6 | isempty(var)
   var = 1;
end
if nargin<7 | isempty(fnms)
   fnms = {'CARS_latest'};
   if isempty(GET_CLIM_def_warn)
      GET_CLIM_def_warn = [0 0];
   end
   if strcmp(prop,'t') | strcmp(prop,'s')
      if GET_CLIM_def_warn(1)==0
	 disp('GET_CLIM uses latest versions of CARS by default')
	 disp('For T and S this is presently cars2009a');
	 GET_CLIM_def_warn(1) = 1;
      end
   elseif GET_CLIM_def_warn(2)==0
      disp('GET_CLIM uses latest versions of CARS by default')
      disp('For nutrients this is presently cars2009');
      GET_CLIM_def_warn(2) = 1;
   end
elseif strcmp(fnms,'best')   
   fnms = {'coast8_06','CARS_latest'};
end
   
if ~iscell(fnms)
   fnms = {fnms};
end

if nargin<4 | isempty(deps)
   [tmp,nc] = clname(prop,fpth,fnms{1});
   deps = nc{'depth'}(:); 
   close(nc);
end

if nargin<8 | isempty(woa); woa = 0; end
if nargin<9 | isempty(fll); fll = 0; end
if nargin<10; fpth = []; end
if nargin<11 | isempty(fns); fns = 2; end

if strcmp(prop,'t')
   nneg = 0;
else
   nneg = 1;
end
if nargin>=12 & ~isempty(config) & isfield(config,'nneg') & ~isempty(config.nneg)
   nneg = config.nneg;
end

vv = interp3_clim(fpth,fnms{1},fll,prop,var,deps,lon,lat,doy,[],fns);
for ii = 2:length(fnms)
   if any(isnan(vv(:)))
      vin = vv;
      vv = interp3_clim(fpth,fnms{ii},fll,prop,var,deps,lon,lat,doy,vin,fns);
   end
end


if any(isnan(vv(:))) & woa>0
   if isempty(strmatch('cars2005',fnms)) & isempty(strmatch('cars2006',fnms))
      kk = find(lat>10 | lon>200 | lon<30);
   else
      kk = find(lat<-70 | lat>26 | (lon>270 & lat>=10));
   end
   if ~isempty(kk)
      disp([num2str(length(kk)) ' profiles used WOA']);
      if isempty(doy)	
	 dok = [];
      else
	 dok = doy(kk);
      end
      
      if strcmp(prop,'si')
	 % WOA prefix "i" used for silicate
	 prop = 'i';
      end
	    
      if woa==1
	 vv(kk,:) = get_woa05_profiles(prop,lon(kk),lat(kk),deps,dok,var)';
      elseif woa==2
	 wdep = dep_csl(deps,1);
	 iw = find(wdep==round(wdep));
	 if length(iw)~=length(deps)
	    nwdep = deps;
	    nwdep(iw) = [];
	    disp('The following depths are not available in WOA98 (is it mapped on');
	    disp(['a smaller set on depth levels: ' num2str(nwdep)]);
	 end
	 if ~isempty(iw)
	    disp([num2str(length(kk)) ' profiles used WOA98']);
	    tmp = get_woa_profiles(prop,lon(kk),lat(kk),wdep(iw),dok);
	    if dims(lon)==1
	       vv(kk,iw) = tmp';
	    else
	       for jj = 1:length(iw)
		  tmp2 = vv(:,:,iw(jj));
		  tmp2(kk) = tmp(jj,:);
		  vv(:,:,iw(jj)) = tmp2;
	       end
	    end
	 end
      end
   end
end  

if nneg
   vv(vv<0) = 0;
end

%------------- End of get_clim -------------------------
