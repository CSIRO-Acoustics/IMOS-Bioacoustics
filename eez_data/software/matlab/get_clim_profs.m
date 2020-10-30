% GET_CLIM_PROFS: Return extracts from CARS, at a set of locations and
%      depths, and optionally day-of-years.   Can also return just the 
%      CARS seasonal anomaly (see 'doy' below).  
%
% INPUT
%  prop    property ('t','s','o','si','n', OR 'p')
%  lon,lat   vectors of NN locations
%  deps    vector of DD depths (in m) [sorted] These need not be CSIRO
%          Standard Levels - ie it will interpolate between CARS levels!
%   OPTIONAL:
%  doy     vector of NN day-of-year corresponding to locations. Note
%          that a seasonal cycle are not provided below a certain depth 
%          (1000m for CARS2000).
%          NOTE: use -doy if want just seasonal component; ie mean=0.
%  fname   CARS version name:  eg 'cars2009a' or 'coast8_06' [def cars_latest]
%  vart    1=seasonal or mean    9=SD   10=RMSMR   11=RMSR  [def 1]
%          (var>=2 disables seasonal evaluation)
% OUTPUT
%  vv      psuedo-casts extracted from CARS, dimensioned [DD,NN]
%  out     index to casts outside region of specified maps
%
% AUTHOR: Jeff Dunn  CSIRO CMAR  Aug 2012
%
% CALLS:  clname  inpolygon  coord2grd  getchunk2
%
% REPLACES get_clim_casts: functionality is unchanged but coded to cope with
%                          changes in Matlab.
%
% RELATED:  'get_clim' - is slower but can get data from multiple climatologies
%           in one go (eg use low res only where no high res), and also can
%           return variables other than the mean or seasonal property fields.
%
% USAGE: [vv,out] = get_clim_profs(prop,lon,lat,deps,doy,fname,vart);

function [vv,out] = get_clim_profs(prop,lon,lat,deps,doy,fname,vart)

[vv,out] = deal([]);
            
if ~isa(lon,'double')
   lon = double(lon);
   lat = double(lat);
end

if nargin<5; doy = []; end
if nargin<6; fname = []; end
if nargin<7 || isempty(vart)
   vart = 1; 
end

lon = lon(:)';
lat = lat(:)';
ncast = length(lat);
if length(doy)==1 && ncast>1
   doy = repmat(doy,[1 ncast]);
else
   doy = doy(:)';
end
nomean = 0;
if ~isempty(doy) & (sum(doy<0)>sum(doy>=0))
   nomean = 1;
   doy = -doy;
end

ndeps = length(deps);
vv = nan(ndeps,ncast);

tcor = -i*2*pi/366;
cpath = [];

cfnm = clname(prop,cpath,fname);
cfnm = [cfnm '.nc'];
gor = ncread(cfnm,'gr_origin');
gsp = ncread(cfnm,'grid_space');
rot = ncread(cfnm,'rotation');
cnrs = ncread(cfnm,'corners');
cldeps = ncread(cfnm,'depth');

[askdeps,icl,jcl] = intersect(deps,cldeps);

if length(icl)<ndeps
   % iint indexes any depths specified which are not in the climatology.
   % ideps are the set of clim depths which straddle each iint depth.
   iint = 1:ndeps;
   iint(icl) = [];
   
   ii = interp1(cldeps,1:length(cldeps),deps(iint));
   ideps = unique([floor(ii); ceil(ii)]);
   ideps(ideps<1 | ideps>length(cldeps)) = [];
   
   askdeps = unique([askdeps(:); cldeps(ideps(:))]);
   
   [tmp,tmp,jicl] = intersect(cldeps(ideps),askdeps);
else
   jicl = [];
end

% Reset jcl so it indexes the std deps in the set of extracted levels 
[tmp,tmp,jcl] = intersect(deps(icl),askdeps);

% If bad location, set to impossible value so it is flagged as out-of-region
lat(isnan(lon) | isnan(lat)) = -999;
out = 1:ncast;
ic = 1:ncast;
if size(cnrs,2)==2
   cnrs = cnrs';
end
if isempty(rot) || rot==0
   if ~isempty(findstr(cfnm,'cars2005')) || ~isempty(findstr(cfnm,'cars2006'))
      out = find(lat<-70 | lat>26 | (lon>270 & lat>=10));
      ic(out) = [];
%   elseif cnrs(1,1)<=-75 && cnrs(2,1)==0 && cnrs(1,3)==90 && cnrs(2,3)==360
%      % global
%      out = [];
   else
      % Apply this test even for global CARS, so impossible locations are detected
      ic = find(lon>=min(cnrs(2,:)) & lon<=max(cnrs(2,:)) & ...
		lat>=min(cnrs(1,:)) & lat<=max(cnrs(1,:)));
      out(ic) = [];
   end
else
   ic = inpolygon(lon,lat,cnrs(2,:),cnrs(1,:));
   ic = find(ic>0);
   out(ic) = [];
end
ic = ic(:)';  

if ~isempty(out)
%   disp(['GET_CLIM_PROFS: ' num2str(length(out)) '/' num2str(ncast) ...
%	 ' locs outside CARS region: ' ...
%	 num2str([range(lat(out)) range(lon(out))])]);
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
     [mn,t2,t3,t4,t5,ix,iy] = getchunk2(prop,askdeps,rng,cpath,fname,-2,vart);
          
     if min(size(mn))<2
        error('GET_CLIM_CASTS - region lacks enough grid points to interpolate')
     end
     X = X+1-min(ix);
     Y = Y+1-min(iy);

     for jj = 1:length(icl)
	vv(icl(jj),ic) = interp2(mn(:,:,jcl(jj)),X,Y,'*linear');
     end
     if ~isempty(jicl)
     	Y = repmat(Y,length(iint),1);
     	X = repmat(X,length(iint),1);
	Z = interp1(askdeps(jicl),1:length(jicl),deps(iint));
     	Z = repmat(Z(:),1,length(ic));
     	vv(iint,ic) = interp3(mn(:,:,jicl),X,Y,Z,'*linear');
     end    
  else
    [mn,an,sa,t4,t5,ix,iy] = getchunk2(prop,askdeps,rng,cpath,fname,2,vart);
     if min(size(mn))<2
	error('GET_CLIM_CASTS - region lacks enough grid points to interpolate')
     end

    if isempty(an)
       %disp('No temporal harmonics available for these depths');
       tdep = 0;
    elseif ndims(an)==2       
       tdep = 1;
    else
       tdep = size(an,3);
    end
    if isempty(sa)
       %disp('No temporal harmonics available for these depths');
       t2dep = 0;
    elseif ndims(sa)==2       
       t2dep = 1;
    else
       t2dep = size(sa,3);
    end
    
    mt = zeros(ndeps,length(ic));
    X = X+1-min(ix);
    Y = Y+1-min(iy);
    if ~nomean
       for jj = 1:length(icl)
	  mt(icl(jj),:) = interp2(mn(:,:,jcl(jj)),X,Y,'*linear');
       end
       if ~isempty(jicl)
	  Y2 = repmat(Y,length(iint),1);
	  X2 = repmat(X,length(iint),1);
	  Z2 = interp1(askdeps(jicl),1:length(jicl),deps(iint));
	  Z2 = repmat(Z2(:),1,length(ic));
	  mt(iint,:) = interp3(mn(:,:,jicl),X2,Y2,Z2,'*linear');
       end    
    end
    if tdep>0
       at = zeros(ndeps,length(ic));
       for jj = 1:length(icl)
	  if jcl(jj)<=tdep
	     at(icl(jj),:) = interp2(an(:,:,jcl(jj)),X,Y,'*linear');
	  end
       end
       if ~isempty(jicl)
	  jicl(jicl>tdep) = [];
       end
       if ~isempty(jicl)	  
	  at(iint,:) = interp3(an(:,:,jicl),X2,Y2,Z2,'*linear');
       end    
       
       if t2dep>0
	  st = zeros(ndeps,length(ic));

	  for jj = 1:length(icl)
	     if jcl(jj)<=t2dep
		st(icl(jj),:) = interp2(sa(:,:,jcl(jj)),X,Y,'*linear');
	     end
	  end
	  if ~isempty(jicl)
	     jicl(jicl>t2dep) = [];
	  end
	  if ~isempty(jicl)	  
	     st(iint,:) = interp3(sa(:,:,jicl),X2,Y2,Z2,'*linear');
	  end    
	  
       end   
    end
    
    % Pre-load vv in case tdep=0 or is less than ndep. 
    vv(:,ic) = mt;
    if tdep>0
       % Replace temporal NaNs with 0 so that these have no effect in computation.
       at(isnan(at)) = 0;

       tdoy = tcor.*repmat(doy(ic),ndeps,1);
       vv(:,ic) = mt + real(at.*exp(tdoy));
       if t2dep>0
	  st(isnan(st)) = 0;
	  vv(:,ic) = vv(:,ic) + real(st.*exp(2*tdoy));
       end
    end
  end
end

%------------- End of get_clim_casts -------------------------
