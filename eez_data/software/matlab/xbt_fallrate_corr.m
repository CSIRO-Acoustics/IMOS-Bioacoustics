% XBT_FALLRATE_CORR  Correct XBT temperature for time-varying fallrate errors
%  as resolved by Wijffels et al 2007. Assumption is that XBTs are already
%  corrected to Hanawa et al 1995 (H95). Note that XBTs outside the
%  correction table time range cannot be corrected. The table was calculated
%  for 1968 to 2005; the algorithm used here allows correction from 1967.5 to
%  2005.5.
%
% INPUT   tz   [ndepth nprof] XBT T(z)
%         z    [ndepth 1] depth grid [same for all profiles]
%      ** typ  [nprof 1]  0=don't correct  1=shallow XBT   2=deep XBT
%         tim  [nprof 1]  decimal days since 1900 of each profile
%         lon  [nprof 1]  longitude (0-360) of each profile
%         lat  [nprof 1]  latitude (north +ve) of each profile
% **If do not know probe type, then can use Wijffels criterion of
%    profile_depth>550m?    ie typ = 1+(pdep>550);  
%
% OUTPUT  tnew - Corrected tz
%
%  6/11/07  Jeff Dunn from Susan Wijffel's algorithm
%
%  Ref:  Changing expendable bathythermograph fall-rates and their impact on
%  estimates of thermosteric sea level rise. S E Wijffels et al 2007 (in prep)
%
% USAGE: tnew = xbt_fallrate_corr(tz,z,typ,tim,lon,lat);

function tnew = xbt_fallrate_corr(tz,z,typ,tim,lon,lat)

tnew = tz;
flipdim = 0;

if size(tz,1)~=length(z)
   if size(tz,2)~=length(z)
      disp('XBT_FALLRATE_CORR:  tz is wrong shape!')
      return
   else
      flipdim = 1;
      tz = tz';      
      tnew = tnew';
   end
end
npro = size(tz,2);

if min(size(z))~=1 || max(size(z))~=size(tz,1)
   disp('XBT_FALLRATE_CORR:  "z" is wrong shape!')
   return
end
z = z(:);

if length(typ)~=npro
   if length(typ)==1
      typ = repmat(typ,[npro 1]);
   else
      disp('XBT_FALLRATE_CORR:  "typ" is wrong shape!')
      return
   end
end
typ = typ(:);

tim = tim(:);   
dectim = 1900 + (tim/365.25);

wmo = ones(size(tim));
wmo(lat<0 & lon<180) = 3;
wmo(lat<0 & lon>=180) = 5;
wmo(lat>0 & lon>180) = 7;

% The corrections are calculated for 2 year bins centred on the given date;
% that is, dcorrs(1,:) are corrections at 1968.0. 
% Load a copy of 
% /home/wijffels/work/global_thermal/slope_corrections_xbt_basin_hanawa.mat

load /home/eez_data/software/data/slope_corrections_xbt_basin_hanawa
dcorrs = squeeze(dcorrs);
dcorrd = squeeze(dcorrd);

% It would seem safe to retrieve an extra year of correction by extending
% the extrapolation range 6 months at either end.
yrgrid = [yrgrid(1)-.5 yrgrid  yrgrid(end)+.5];
dcorrs = [dcorrs(1,:); dcorrs; dcorrs(end,:)];
dcorrd = [dcorrd(1,:); dcorrd; dcorrd(end,:)];

ncorr = 0;

for ptyp = 1:2
   ii = find(dectim>=yrgrid(1) & dectim<=yrgrid(end) & typ==ptyp);
   ncorr = ncorr+length(ii);
   if ptyp==1
      dcorr = dcorrs;
   else
      dcorr = dcorrd;
   end   
   
   for basin = 1:4
      ip = ii(wmo(ii)==basin);
      dze = z*interp1q(yrgrid',squeeze(dcorr(:,basin)),dectim(ip))';      
      for jj = 1:length(ip)
	 % xbt's are reporting depths that are too deep - need to reduce.
	 tnew(:,ip(jj)) = interp1q(z-dze(:,jj),squeeze(tz(:,ip(jj))),z); 
      end
   end
end

disp([num2str(ncorr) ' out of ' num2str(npro) ' XBTs corrected']);

if flipdim
   tnew = tnew';
end

%----------------------------------------------------------------------------
