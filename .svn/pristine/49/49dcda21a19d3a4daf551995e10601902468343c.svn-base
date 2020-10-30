% ETOPO2V2: Extract heights from etopo2v2 datafile, for given lats & lons.
%           This dataset is global, 2-minute resolution, cell centred.
%
% INPUT:
%    lat - matrix of lats, -90-90deg
%    lon - matrix of lons,   0-360deg
%
% OUTPUT:
%    z - in metres, +ve upwards (so heights, not depths). 
%
% NOTE the 2-min dataset is just sampled rather than using 2D interpolation.
%
% USAGE: z = etopo2v2(lats,lons)

function z = etopo2v2(lats,lons)

fname = path_pc_or_nix('netcdf-data/ETOPO2v2c_f4_cmar.nc');

z = repmat(NaN,size(lats));
lons = lons(:);
lats = lats(:);

if any(lons<0)
   if any(lons>180) || any(lons<-180)      
      disp('ETOPO2: Longitudes wrong, use 0-360 or -180 to 180!')
      return      
   else
      % Do nothing - apparently longitudes already -180 to 180. 
   end
elseif any(lons>360)
   disp('ETOPO2: Longitudes wrong, use 0-360 or -180 to 180!')
   return
else
   % Expect lons to be 0-360, so convert to -180-180 to suit input file
   kk = find(lons>180);
   lons(kk) = lons(kk)-360;
end
if any(abs(lats)>89.983)
   if any(abs(lats)>90)
      disp('ETOPO2V2: Latitude outside range -90 to 90');
      return
   else
      ii = find(abs(lats)>89.983 & abs(lats)<=90);
      lats(ii) = sign(lats(ii))*89.983;
   end
end

ix = 1 + round((lons+179.983)*30);
iy = 1 + round((lats+89.983)*30);
if any(ix==0)
   ix(ix==0) = 1;
end
if any(ix==10801)
   ix(ix==10801) = 10800;
end

jx = [min(ix) max(ix)];
jy = [min(iy) max(iy)];

for kx = jx(1):500:jx(2)
   mm = find(ix>=kx & ix<(kx+500) & ~isnan(iy));
   if ~isempty(mm)   
      hh = getnc(fname,'height',[jy(1) kx],[jy(2) min(kx+499,10800)]);
      if min(size(hh))==1
	 % Getnc stupidly returns an [N 1] vector, irrespective of whether
	 % asked for a [1 N] or an [N 1] slice, so have to treat each case
	 % separately.
	 if max(1+iy(mm)-jy(1))>1
	    ii = 1+iy(mm)-jy(1);
	 else
	    ii = 1+ix(mm)-kx;
	 end
      else
	 ii = sub2ind(size(hh),1+iy(mm)-jy(1),1+ix(mm)-kx);	 
      end
      z(mm) = hh(ii);
   end
end
      

% ------------ End of etopo2v2.m -------------------
