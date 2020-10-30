% ATDAYPOS  Given lats,lons and times, return values according to grid and 
%           temporal functions.
% INPUT:  
%   lats, lons, doy:  vectors of required positions and day-of-year
%   xgrid,ygrid, mn:  [ny,nx] long, lat and mean map
%   an,sa          :  [ny,nx] optional annual and semi-annual harmonics
%   cap [optional] : minimum value. Best fit seasonal harmonics may sometimes 
%              create -ve values for some part of the year (which for most
%              properties is not physically reasonable.) Eg for nutrients
%              could use 0, for MLD might use 3?, for T use -3.5?
%
% Assumes: Uniformly spaced grid, with element 1,1 being in the SW corner. 
%
% JRD 23/4/96         Change to faster V7 code, 26/2/08
%
% USAGE: vals = atdaypos(lats,lons,doy,xgrid,ygrid,mn,an,sa,cap)

function vals = atdaypos(lats,lons,doy,xgrid,ygrid,mn,an,sa,cap)

if nargin<9
   cap = [];
end

t1 = 0;
gotan = nargin > 6 && ~isempty(an);
if gotan
   t1 = any(isnan(an(:)));
end
t2 = 0;
gotsa = nargin > 7 && ~isempty(sa);
if gotsa
  t2 = any(isnan(sa(:)));
end

if min(size(xgrid))==1 & min(size(mn)) > 1
   [xgrid,ygrid] = meshgrid(xgrid,ygrid);
end

gotnan = any(isnan(lats) | isnan(lons) | isnan(doy));
if gotnan
   [niy,nix] = size(lats);
   jgd = find(~isnan(lats) & ~isnan(lons) & ~isnan(doy));
   lats = lats(jgd);
   lons = lons(jgd);
   doy = doy(jgd);
end

% If there are no NaNs in the harmonics we can get away with the much simpler
% and more accurate method below:
    
if ~(t1 || t2)

  vals = interp2(xgrid,ygrid,mn,lons,lats,'*bilinear');

  if gotan
    ann = interp2(xgrid,ygrid,an,lons,lats,'*bilinear');
    vals = vals + real(ann.*exp(-i*doy*2*pi/366));
  end
  if gotsa
    saa = interp2(xgrid,ygrid,sa,lons,lats,'*bilinear');
    vals = vals + real(saa.*exp(-i*doy*4*pi/366));
  end

else
   if max(size(doy))==1 & max(size(lats))>1
      doy = repmat(doy,size(lats));
   end
   
   an(isnan(an)) = 0;
   if t2
      sa(isnan(sa)) = 0;
   end
   
   % For pre-version7 Matlab, use:
   %ii = find(isnan(an));
   %an(ii) = zeros(size(ii));
   %if t2
   %   ii = find(isnan(sa));
   %   sa(ii) = zeros(size(ii));
   %end
  
  [mg,ng] = size(xgrid);
  vals = NaN*ones(length(lats),1);

  vbase = interp2(xgrid,ygrid,mn,lons,lats,'*bilinear');

  gsp = abs(xgrid(1,1)-xgrid(1,2));
  orgn = [xgrid(1,1)-gsp/2 ygrid(1,1)-gsp/2];

  if ~gotan
    vals = vbase;
  elseif ~gotsa
    for ii = 1:length(lats)
      elr = 1+floor((lons(ii)-orgn(1))/gsp);
      elc = 1+floor((lats(ii)-orgn(2))/gsp);
      if( elr>=1 & elr<=ng & elc>=1 & elc<=mg)
	vals(ii) = vbase(ii) + real(an(elc,elr)*exp(-i*doy(ii)*2*pi/366));
      end
    end
  elseif gotsa
    for ii = 1:length(lats)
      elr = 1+floor((lons(ii)-orgn(1))/gsp);
      elc = 1+floor((lats(ii)-orgn(2))/gsp);
      if( elr>=1 & elr<=ng & elc>=1 & elc<=mg)
	vals(ii) = vbase(ii) ...
	    + real(an(elc,elr)*exp(-i*2*pi/366*doy(ii))) ...
	    + real(sa(elc,elr)*exp(-i*4*pi/366*doy(ii)));
      end
    end
  end
  
end

if gotnan
   tmp = vals;
   vals = repmat(nan,[niy nix]);
   vals(jgd) = tmp;
end

if ~isempty(cap)
   vals(vals<cap) = cap;
end

% -------------- End of atdaypos ----------------------

