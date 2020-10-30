% GET_CLIM_CASTS: Return extracts from CARS, at a set of locations and
%      depths, and optionally day-of-years.   Can also return just the 
%      CARS seasonal anomaly (see 'doy' below).  SUPERCEDES getsection
%
% INPUT
%  prop    property ('t','s','o','si','n', OR 'p')
%  lon,lat   vectors of NN locations
%  deps    vector of DD depths (in m) [sorted but not nec. contiguous]
%   OPTIONAL:
%  doy     vector of NN day-of-year corresponding to locations. Note
%          that a seasonal cycle are not provided below a certain depth 
%          (1000m for CARS2000).
%          NOTE: use -doy if want just seasonal component; ie mean=0.
%  fname   CARS version name:  eg 'cars2006a' or 'coast8_06' [def cars2006]
%  woa     1=WOA05  2=WOA98   used outside of CARS region  [def 0]
%  fll     1=values at all depths outside of land areas (with selected
%          files only) [default 0]
%  vart    1=seasonal or mean  9=SD  10=RMSMR  11=RMSR  [def 1]
%          (var>=2 disables seasonal evaluation)
% OUTPUT
%  vv      psuedo-casts extracted from CARS, dimensioned [DD,NN]
%  out     index to casts outside region of specified maps
%
% AUTHOR: Jeff Dunn  CSIRO DMR  May 1998
% $Id: get_clim_casts.m,v 1.6 2001/06/07 00:00:59 dun216 Exp dun216 $
%
% CALLS:  clname  inpolygon  coord2grd  getchunk  get_woa_profiles
% RELATED:  'get_clim' - is slower but can get data from multiple climatologies
%           in one go (eg use low res only where no high res), and also can
%           return variables other than the mean or seasonal property fields.
%
% USAGE: [vv,out] = get_clim_casts(prop,lon,lat,deps,doy,fname,woa,fll,vart);

function [vv,out] = get_clim_casts(prop,lon,lat,deps,doy,fname,woa,fll,vart);

% MODS  30/5/08  Add access to WOD05

[vv,out] = deal([]);

ncquiet;

if ~isa(lon,'double')
   lon = double(lon);
   lat = double(lat);
end

if nargin<5; doy = []; end
if nargin<6; fname = []; end
if nargin<7 | isempty(woa); woa = 0; end
if nargin<8 | isempty(fll); fll = 0; end
if nargin<9 | isempty(vart)
   vart = 1; 
   doy = [];
end

ndep = length(deps);

lon = lon(:)';
lat = lat(:)';
ncast = length(lat);
if length(doy)==1 & ncast>1
   doy = repmat(doy,[1 ncast]);
else
   doy = doy(:)';
end
nomean = 0;
if ~isempty(doy) & (sum(doy<0)>sum(doy>=0))
   nomean = 1;
   doy = -doy;
end

vv = repmat(nan,ndep,ncast);

tcor = -i*2*pi/366;
cpath = [];

cfnm = clname(prop,cpath,fname);
ncf = netcdf([cfnm '.nc'],'nowrite');
gor = ncf{'gr_origin'}(:);
gsp = ncf{'grid_space'}(:);
rot = ncf{'rotation'}(:);
cnrs = ncf{'corners'}(:);
close(ncf);


out = 1:ncast;
ic = 1:ncast;
if isempty(rot) | rot==0
   if isempty(findstr(cfnm,'cars2005')) & isempty(findstr(cfnm,'cars2006'))
      ic = find(lon>=min(cnrs(2,:)) & lon<=max(cnrs(2,:)) & ...
		lat>=min(cnrs(1,:)) & lat<=max(cnrs(1,:)));
      out(ic) = [];
   else
      out = find(lat<-70 | lat>26 | (lon>270 & lat>=10));
      ic(out) = [];
   end
else
   ic = inpolygon(lon,lat,cnrs(2,:),cnrs(1,:));
   ic = find(ic>0);
   out(ic) = [];
end
ic = ic(:)';  

if ~isempty(out)
   disp(['GET_CLIM_CASTS: ' num2str(length(out)) '/' num2str(ncast) ...
	 ' locs outside CARS region: ' ...
	 num2str([range(lat(out)) range(lon(out))])]);
end


if ~isempty(ic) 
  % Auto-set an efficient range, but guard against degenerate ones which 
  % would not provide enough grid points for interpolation.

  rng = [floor((min(lon(ic))-.1)/gsp(2))*gsp(2) ...
	  ceil((max(lon(ic))+.1)/gsp(2))*gsp(2) ...
	 floor((min(lat(ic))-.1)/gsp(1))*gsp(1) ...
	  ceil((max(lat(ic))+.1)/gsp(1))*gsp(1)];
  
  % Convert position to index-coords into the climatology chunk, so can
  % use abbreviated form for interp3 (ie not supply meshgrid). Because deeper
  % NaNs wipe out estimates up to and including at the layer above, we shift Z
  % to be just above the layer so they are uneffected by NaNs below. (However, 
  % this would then put layer 1 outside of the grid and hence lose it, so we
  % 2D-interpolate it separately!)

  [X,Y] = coord2grd(lon(ic),lat(ic),gor(2),gor(1),gsp(2),gsp(1),rot);

  if isempty(doy)
     [mn,t2,t3,t4,t5,ix,iy] = getchunk(prop,deps,rng,cpath,fname,-2,fll,vart);
     if isempty(mn)
	return
     end
     
     X = X+1-min(ix);
     Y = Y+1-min(iy);
     if min(size(mn))<2
	error('GET_CLIM_CASTS - region lacks enough grid points to interpolate')
     end
     for ii = 1:ndep
	vv(ii,ic) = interp2(mn(:,:,ii),X,Y,'*linear');
     end
     %vv(1,ic) = interp2(mn(:,:,1),X,Y,'*linear');
     %if ndep>1
     %	Y = repmat(Y,ndep-1,1);
     %	X = repmat(X,ndep-1,1);
     %	Z = repmat((2:ndep)'-.0001,1,length(ic));
     %	vv(2:ndep,ic) = interp3(mn,X,Y,Z,'*linear');
     %end    
  else
    [mn,an,sa,t4,t5,ix,iy] = getchunk(prop,deps,rng,cpath,fname,2,fll,vart);
     if isempty(mn)
	return
     end

    if isempty(an)
       %disp('No temporal harmonics available for these depths');
       tdep = 0;
    elseif ndims(an)==2       
       tdep = 1;
    else
       tdep = ndep;
    end
    if isempty(sa)
       semian = 0;
    else
       semian = 1;
    end
    
    X = X+1-min(ix);
    Y = Y+1-min(iy);
    if nomean
       mt = zeros(ndep,length(ic));
    else
       mt = interp2(mn(:,:,1),X,Y,'*linear');
    end
    if tdep>0
       at = interp2(an(:,:,1),X,Y,'*linear');
       if semian; st = interp2(sa(:,:,1),X,Y,'*linear'); end
    end
    
    if ndep>1
      Y = repmat(Y,ndep-1,1);
      X = repmat(X,ndep-1,1);
      Z = repmat((2:ndep)'-.0001,1,length(ic));
      if ~nomean
	 mt = [mt; interp3(mn,X,Y,Z,'*linear')];
      end
      clear mn
      if tdep>1
	 at = [at; interp3(an,X,Y,Z,'*linear')];
	 clear an
	 if semian; st = [st; interp3(sa,X,Y,Z,'*linear')]; end
      end
    end

    % Pre-load vv in case tdep=0 or is less than ndep. 
    vv(:,ic) = mt;
    if tdep>0
       % Replace temporal NaNs with 0 so that these have no effect in computation.
       kk = find(isnan(at));
       if ~isempty(kk), at(kk) = zeros(size(kk)); end
       if semian
	  kk = find(isnan(st));
	  if ~isempty(kk), st(kk) = zeros(size(kk)); end
       end
       tdoy = tcor.*repmat(doy(ic),tdep,1);

       if semian
	  vv(1:tdep,ic) = mt(1:tdep,:) + real(at.*exp(tdoy)) + real(st.*exp(2.*tdoy));
       else
	  vv(1:tdep,ic) = mt(1:tdep,:) + real(at.*exp(tdoy));
       end
    end
  end
end


if woa>0 & ~isempty(out)
   if isempty(doy)	
      dok = [];
   else
      dok = doy(out);
   end      
   if strcmp(prop,'si')
      % WOA prefix "i" used for silicate
      prop = 'i';
   end
   if woa==1
      vv(:,out) = get_woa05_profiles(prop,lon(out),lat(out),deps,dok,vart);	   
   elseif woa==2
      wdep = dep_csl(deps,1);
      iw = find(wdep==round(wdep));
      if length(iw)~=ndep
	 nwdep = deps;
	 nwdep(iw) = [];
	 disp('The following depths are not available in WOA98 (is it mapped on');
	 disp(['a smaller set on depth levels: ' num2str(nwdep)]);
      end
      if ~isempty(iw)
	 vv(iw,out) = get_woa_profiles(prop,lon(out),lat(out),wdep(iw),dok,vart);
      end
   end
   out = [];      
end   

%------------- End of get_clim_casts -------------------------
