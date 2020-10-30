% XBT_FALLRATE_CORR_13  Correct XBT temperature for time-varying fallrate errors
%  as resolved by Cowley et al 2013. Assumption is that XBTs are already
%  corrected to Hanawa et al 1995 (H95). Note that XBTs outside the
%  correction table time range cannot be corrected. The table was calculated
%  for 1967.5 to 2010.5
%
% INPUT   tz   [ndepth nprof] XBT T(z)
%         z    [ndepth 1] depth grid [same for all profiles]
%         deep  [nprof 1]  0=don't correct  1=shallow XBT   2=deep XBT  
%              (See /home/dunn/eez_data/csiro_therm_archive/iota_all/new_probe_clasif.m)
%         tim  [nprof 1]  decimal days since 1900 of each profile
%         TSK  [nprof 1]  1=TSK probe 
%         Tnum  [nprof 1] eg 4 => T4 probe  (presently only used to prevent
%                                            T5 correction)
%
% OUTPUT  tnew - Corrected tz
%         znew - new Z values, if NOT interpolating back onto original z
%
% NOTE  tz will be interpolated back to original z if only one output
%       argument is used. This particularly suits standard depth data.
%
%  18/4/2013 Jeff Dunn
%
%  Ref:  http://www.nodc.noaa.gov/OC5/XBT_BIAS/cowley.html
%    Application of corrections to original XBT data
%    based on Cowley, R., S. Wijffels, L. Cheng, T. Boyer, S. Kizu: 
%    Biases in Expendable BathyThermograph data: a new view based on
%    historical side-by-side comparisons, 
%    accepted by the Journal of Atmospheric and Oceanic Technology. 
%
% USAGE: [tnew,znew] = xbt_fallrate_corr(tz,z,deep,tim,TSK,Tnum);

function [tnew,znew] = xbt_fallrate_corr_13(tz,z,deep,tim,TSK,Tnum)

origz = (nargout==1);

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

if origz
   znew = [];
else
   znew = repmat(z,[1 npro]);
end

if length(deep)~=npro
   if length(deep)==1
      deep = repmat(deep,[npro 1]);
   else
      disp('XBT_FALLRATE_CORR:  "deep" is wrong shape!')
      return
   end
end
deep = deep(:);

if nargin<5 || isempty(TSK)
   % If TSK / Sippican distinction not made then assume Sippican
   TSK = zeros(size(deep));
end

if nargin<6 || isempty(Tnum)
   % If no Tnum then just use "deep"
   Tnum = [];
else
   deep(Tnum==5) = 0;         % No correction available for T5
end

tim = tim(:);   
dectim = 1900 + (tim/365.25);

% Correction group:  1 = T4/T6,  2 = T7/DB,  3 = TSK T6,  4 = TSK T7
ncor = zeros(npro,1);
ncor(deep==1 & ~TSK) = 1;
ncor(deep==2 & ~TSK) = 2;
ncor(deep==1 &  TSK) = 3;
ncor(deep==2 &  TSK) = 4;

% Load: alpha beta deltaT yr
load /home/eez_data/software/data/xbt_corrections_cowley13

ncorr = 0;

for ctyp = 1:4
   ii = find(dectim>=yr(1) & dectim<=yr(end) & ncor==ctyp);
   ncorr = ncorr+length(ii);

   A = interp1q(yr,alpha(:,ctyp),dectim(ii));
   B = interp1q(yr,beta(:,ctyp),dectim(ii));
   D = interp1q(yr,deltaT(:,ctyp),dectim(ii));

   for jj = 1:length(ii)
      zcor = z*(1-A(jj))-B(jj);
      tzd = tz(:,ii(jj))-D(jj);

      if origz
	 t = interp1q(zcor,tzd,z);
      
	 % Shifting the bottom value up means that value cannot be interpolated
	 % back to original Z, so we have to find and extrapolate those cases.
	 if any(isnan(t) & ~isnan(tzd))
	    kk = find(isnan(t) & ~isnan(tzd));
	    kk(kk<3) = [];	 
	    if any(kk)
	       t(kk) = tzd(kk) + (tzd(kk)-tzd(kk-1)).*(z(kk)-zcor(kk))./(zcor(kk)-zcor(kk-1));
	    end
	 end
	 
	 tnew(:,ii(jj)) = t;
      else
	 tnew(:,ii(jj)) = tzd;
	 znew(:,ii(jj)) = zcor;
      end	 
   end
end

disp([num2str(ncorr) ' out of ' num2str(npro) ' XBTs corrected']);

if flipdim
   tnew = tnew';
   znew = znew';
end

%----------------------------------------------------------------------------
