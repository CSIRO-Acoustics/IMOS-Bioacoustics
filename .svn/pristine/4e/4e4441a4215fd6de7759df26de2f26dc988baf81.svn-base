% GA9SECBATH  Returns either full-resolution chunk or values at locations, 
%  from GA 9-sec bathymetry. In the latter case, returns nearest values to 
%  locations, rather than interpolating.
%
% INPUTS
%      x,y     locations at which data required
%            OR
%      region   either [w e s n]     OR     [x1 y1; x2 y2; x3 y3; ... xn yn]
%
%   OUTPUTS
%      deps    [same shape as inputs] depths (m), -ve upwards (+ve depth).
%      x,y     IF using "region", locations of full res values within region.	
%      vers    1= 2005   2=2009   [default 2]
%
% Warning: The region requested should be kept very small to avoid 
%          crashing Matlab.
% NOTE the 9-sec dataset is just sampled rather than using 2D interpolation.
%
% Jeff Dunn CMAR 13/8/07
% MODS:  5/3/2010  Software unchanged, but switched to a different build of
%        the ga_9sec_2009.nc    JRD
%
% EXAMPLES:    deps = ga9secbath(x,y);
%      OR:     [deps,x,y] = ga9secbath([],[],[152 154 -34 -28]);
%      OR:     px = [152 154 155 154 153]; py = [-34 -34 -32 -30 -30];
%              [deps,x,y] = ga9secbath([],[],[px(:) py(:)]);
%
% USAGE: [deps,x,y] = ga9secbath(x,y,region,vers);

function [deps,x,y] = ga9secbath(xin,yin,region,vers);

deps = []; x = []; y = [];

if nargin<4 || isempty(vers)
   vers = 2;
end

padd = ~exist('path_pc_or_nix','file');
if padd
   addpath /home/eez_data/software/matlab/
end

if vers==1
   % ga_9sec_bath is on a grid of X = 92:.0025:171.9975, Y=-59.9975:.0025:-8
   fname = path_pc_or_nix('netcdf-data/ga_9sec_bath.nc');
   x1 = 92;  x2 = 179.9975;
   nx = 35200;
   y1 = -59.9975;  y2 = -8;
else
   % ga_9sec_2009 is X = 91.99875:.0025:171.99625, Y=-60.00125:.0025:-8.00375
   fname = path_pc_or_nix('datalib/bathymetry/GA_BATHYTOPO_200906/data/netcdf/bath_ga_2009.nc');
   x1 = 91.99875;   x2 = 171.99625;
   nx = 32000;
   y1 = -60.00125;  y2 = -8.00375;
end

if padd
   rmpath /home/eez_data/software/matlab/
end

if nargin==2
   region = [];
elseif nargin<2
   disp('ga9secbath: At least 2 inputs required')
   return
elseif nargin>=3 && ~isempty(yin) && ~isempty(region)    
   disp('ga9secbath: Either use inputs 1&2, OR input 3')
   return
end

if isempty(region)
   deps = repmat(NaN,size(xin));
   xin = xin(:);
   yin = yin(:);

   if any(yin<y1 | yin>y2 | xin<x1 | xin>x2)
      ii = find(yin>=y1 & yin<=y2 & xin>=x1 & xin<=x2);
   else
      ii = 1:length(xin);
   end
   ix = 1 + round((xin(ii)-x1)*400);
   iy = 1 + round((yin(ii)-y1)*400);

   jx = [min(ix) max(ix)];
   jy = [min(iy) max(iy)];
   
   for kx = jx(1):500:jx(2)
      mm = find(ix>=kx & ix<(kx+500));
      if ~isempty(mm)   
	 hh = getnc(fname,'height',[jy(1) kx],[jy(2) min(kx+499,nx)]);
	 jj = sub2ind(size(hh),1+iy(mm)-jy(1),1+ix(mm)-kx);
	 deps(ii(mm)) = hh(jj);
      end
   end

else
   X = x1:.0025:x2;
   Y = y1:.0025:y2;
   px = [];
   if ~all(size(region)==[1 4])
      if size(region,2)==2
	 px = region(:,1);
	 py = region(:,2);	 
	 region = [min(px) max(px) min(py) max(py)];
      else
	 disp('ga9secbath:  "region" wrong shape.')
	 return
      end
   end
   
   ix = 1 + round((region([1 2])-x1)*400);
   iy = 1 + round((region([3 4])-y1)*400);
   if ix(1)<1; ix(1) = 1; end
   if ix(2)>nx; ix(2) = nx; end
   if iy(1)<1; iy(1) = 1; end
   if iy(2)>20800; iy(2) = 20800; end

   deps = getnc(fname,'height',[iy(1) ix(1)],[iy(2) ix(2)]);
   [x,y] = meshgrid(X(ix(1):ix(2)),Y(iy(1):iy(2)));

   if ~isempty(px)
      ip = inpolygon(x,y,px,py);
      x = x(ip);
      y = y(ip);
      deps = deps(ip);
   end   
end
      
deps = -deps;

% ------------ End of ga9secbath.m -------------------
