% GET_ALT_XY: Return altimeter data (interpolated in space and time from
%  AVISO 1/4 deg, 10 day, gridded dataset) at given locations and times. 
%
% SEE ALSO   get_alt_xy_nrt.m  for access to locally mapped fields
%
% TEMPORAL COVERAGE:  Oct 1992 - recent (periodically updated) 
% 
% INPUT
%  x,y  - vector or matrices of locations
%  tim  - Unix time of required data; eg greg2time([2001 12 31 0 0 0]).
%         Either single or one for every location.
%  thr  - [optional] Maximum % error threshold - applied to the "quadratic
%         mapping error". If omitted, no restriction applied. (Applied after
%         also interpolating error fields to the x-y locations, rather than
%         using this to screen the original gridded altimeter fields prior
%         to interp "alt" values.)
%  opt  - Dataset to use.  [default = 6]   
%       1  sat_alt_msla_global.nc           0.00 -> 359.00   -81.91 ->  81.97
%       2  sat_alt_msla_aus.nc             80.00 -> 180.00   -82.00 ->  14.74
%       3  sat_alt_msla_southern_ocean.nc   0.00 -> 359.67   -82.00 -> -30.25
%       4  sat_alt_msla_south_indian.nc     0.00 -> 165.00   -59.98 -> -20.01
%       5  sat_alt_msla_south_pacific.nc  125.00 -> 300.00   -82.00 ->  14.74
%       6  sat_alt_msla_global_full.nc      0.00 -> 359.67   -82.00 ->  81.97
%
%      As at Feb 2010 they span 14/10/92 - 22/07/09.
%      Completely new versions released 25/11/2009, including name change
%      from "tp_ers" to "sat_alt" to reflect the fact that data from several
%      satellites is included. The combination varies with time, an 
%      indication of what data can be included is:
%                TOPEX/Poseidon: 1992-2005
%                ERS-1:    1992-1993 & 1995-1996
%                ERS-2:    1995-2003
%                GFO:      2000-2008
%                Jason-1:  2002-present
%                Envisat:  2003-present
%                Jason-2:  2008-present
%       sat_alt_msla_global.nc is on a 1 degree grid, all the rest are on
%       a 1/3rd degree Mercator grid.
%
% OUTPUT
%  alt  - alt at x,y,t locations (nan where no data available)
%  dset - no longer used
%  err  - raw error fields (as a %) (max of the 2 straddling timeslices)
%
% Jeff Dunn  CSIRO CMR 5/4/02  
%
% MODS:  See below
%
% USAGE: [alt,dset,aerr] = get_alt_xy(x,y,tim,thr,opt);

function [alt,dset,aerr] = get_alt_xy(x,y,tim,thr,opt)

% MODS: 
%  24/3/04  Allow for datasets 2-4
%  22/2/06 New files do not have variable Alt_flag, so reading it is
%        now disabled. 
%  24/4/08  Add dataset 7
%  25/2/10  All new files, now in /home/datalib/

if nargin<5 | isempty(opt)
   opt = 6;
end

if nargin<4
   thr = [];
end
geterr = (nargout>2 | ~isempty(thr));

if max(size(tim))==1
   tim = repmat(tim,size(x));
end

% Create output variables and set up loop control
alt = repmat(nan,size(x));
aerr = repmat(nan,size(x));
dset = [];

apth = path_pc_or_nix('datalib/observations/rs/altimetry/gridded/MSLA/');
switch opt
  case 1
    %  Global 1 degree: 1 deg lon; lat 1 deg at Eq to ~.14 deg at ~82S/N.
    % 7 day interval
    infl = [apth 'sat_alt_msla_global'];
    tim0 = greg2time([1992 10 14 0 0 0]);
    ming = 1.1;
    nda = 7;
  case 2 
    %  Aus: 1/3 deg 80 - 180E,  lat 1/3 deg at ~15S, ~.05 deg at 82S,
    %  range 82S to 14.74N  7 day interval
    infl = [apth 'sat_alt_msla_aus'];
    tim0 = greg2time([1992 10 14 0 0 0]);
    ming = .34;
    nda = 7;
  case 3
    %  Southern Ocean: 1/3 degree 0-360E, lat .05 @ 82S to .29 @ 30.25S 
    infl = [apth 'sat_alt_msla_southern_ocean'];
    tim0 = greg2time([1992 10 14 0 0 0]);
    ming = .34;
    nda = 7;
  case 4
    %  South Indian: 1/3 deg 0-165E,  lat .31 @ 20S to .17 @ 60S.
    infl = [apth 'sat_alt_msla_south_indian'];
    tim0 = greg2time([1992 10 14 0 0 0]);
    ming = .34;
    nda = 7;
  case 5
    %  South Pacific: 125-300E, 82S-14.74N, lat .33 @ Eq - .05 @ 82S.
    %  7 day interval
    infl = [apth 'sat_alt_msla_south_pacific'];
    tim0 = greg2time([1992 10 14 0 0 0]);
    ming = .34;
    nda = 7;
  case 6
    %  ~1/3 deg Global: 0-360E, 82S-82, lat .33 @ Eq - .05 @ 82S.
    %  7 day interval
    infl = [apth 'sat_alt_msla_global_full'];
    tim0 = greg2time([1992 10 14 0 0 0]);
    ming = .34;
    nda = 7;
  otherwise
    error(['Do not know option ' num2str(opt)]);
end

% 22/2/06 This missing from some files, so disable read to prevent crash
%  aflg = getnc(infl,'Alt_flag');

lo = getnc(infl,'lon');
la = getnc(infl,'lat');
atim = getnc(infl,'time');
atim = atim+tim0;
ntim = length(atim);
atim(end) = atim(end)+.0001;   %little fudge in case requested date falls on
                               %the last day of data. 

for ii = 1:(ntim-1)
   jj = find(tim>=atim(ii) & tim<atim(ii+1));

   if ~isempty(jj)
      % Extract a minimum rectangle and interpolate height and error values
      ix = find(lo>(min(x(jj))-ming) & lo<max(x(jj))+ming);
      iy = find(la>(min(y(jj))-ming) & la<max(y(jj))+ming);

      if length(ix)>1 & length(iy)>1
	 % then we must be inside the region covered, and so can do a 3D interp.
	 hgt = getnc(infl,'height',[ii iy(1) ix(1)],[ii+1 iy(end) ix(end)]);
	 tjj = tim(jj)-atim(ii);
	 alt(jj) = interp3(lo(ix),la(iy),[0 nda],shiftdim(hgt,1),x(jj),y(jj),tjj);

	 if geterr
	    err = getnc(infl,'mapping_error',[ii iy(1) ix(1)],[ii+1 iy(end) ix(end)]);
	    err = squeeze(max(shiftdim(err,1),[],3));
	    aerr(jj) = interp2(lo(ix),la(iy),err,x(jj),y(jj));
	 end
      end
   end
end
   
% If screening, nan results where interpolated errors to high
% We screen after interpolation because if we screened initially to put nans 
% in the gridded fields, then each nan infects all interpolation between it
% and its 6 neighbours - which is pretty severe!

if ~isempty(thr)
   jj = find(aerr>thr);
   alt(jj) = repmat(nan,size(jj));
end

%---------------------------------------------------------------------------
