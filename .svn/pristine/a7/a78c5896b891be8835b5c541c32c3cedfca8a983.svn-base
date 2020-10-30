% S_FROM_DIRECT_TS  Given T casts, return matching S casts from direct S(T)
%       climatology. 
%
%  NOTE:  T above range of S(T) climatology (32C) is capped at 32 to obtain
%         a value.
%
%  NOTE:  In some regions the relationship of S to T can be weak or variable
%         (in all or a part of the T domain.)  Some measure of this is seen
%         in the output arguments "rmsr" and "rmsmr", which approximate a
%         standard deviation of S(T) at each (T,x,y) point. 
%
% INPUTS
%  lo    longitude of each profile
%  la    latitude    "       "
%  doy   day-of-year of each profile (ie 1-366), 
%        OR doy=[] to disable seasonal calc
%  t     temperature profiles   [ndep,ncast]
%  vers  [optional] version number of TS climatology to use
%          3 = May 2007 (with Argo) 
%         10 = latest version (presently same as last listed, above) [Default]
%  fnm   [optional]  Can use INSTEAD of specifying "vers" if don't want a 
%        standard climatology. Full-path name of the alternate netCDF file 
%        (but leave off '.nc') 
%
% OUTPUTS
%   NOTE - if any locations outside the climatology domain then -
%              a) "s" will NOT be the same size as "t"
%              b) "outofreg" will be non-empty
%
%  s    [ndep,NC2]  NC2 is number of casts within region of climatology 
%  outofreg   index to out-of-region casts
%  rmsr       same size as s, the RMS of residuals of S wrt seasonal S(T)
%  rmsmr      same size as s, the RMS of residuals of S wrt mean S(T), hence
%             includes some measure of the seasonal vairability.
%
% Author: Jeff Dunn 4 Nov 2005
%
% Documentation: www.marine.csiro.au/~dunn/ts_clim/index.html
% Climatology created in ~dunn/CARS/p_vs_p/
% Software now in eez_data/software/loess_mapping/p_vs_p/
% Climatology stored in eez_data/ts_clim/
%
% USAGE: [s,outofreg,rmsr,rmsmr] =  s_from_direct_ts(lo,la,doy,t,vers,fnm);


% WARNING:  This algorithm looks tacky and unnecessary. It is also about 150x
%  faster than the obvious interp3 approach, so stick with it!
%
% Mods:  17/5/07  added rr1 rr2 retrieval option
%        29/8/07  corrected to cope with NaNs being in different places in
%                 mean and harmonics.
%        21/5/15  Tweaking to bring code into BOA suite

function [ss,outofreg,rr1,rr2] =  s_from_direct_ts(lo,la,doy,tt,vers,fnm)


[ss,outofreg,rr1,rr2] = deal([]);


if nargin<6 || isempty(fnm)   
   if nargin<5 || isempty(vers)
      vers = 10;
   end
else
   if nargin>=5 && ~isempty(vers)
      disp('S_FROM_DIRECT_TS: Both vers AND fnm specified - ignoring vers!')  
   end
   vers = 0;
end

val_vers = [0 3 10];
if ~any(val_vers==vers)
   warning('S_FROM_DIRECT_TS: No such version! Using default!')
   vers = 10;
end

if vers==3
   fnm = path_pc_or_nix('eez_data/ts_clim/ts_clim07');
else
   fnm = path_pc_or_nix('eez_data/ts_clim/ts_clim');
end

lo = lo(:)';
la = la(:)';
if ~isempty(doy)
   doy = doy(:)';
   ann = 1;
else
   ann = 0;
end
 
if size(tt,2)~=length(lo)
   if size(tt,1)==length(lo)
      tt = tt';
   else
      error('S_FROM_DIRECT_TS: wrong dimensions for "t"')
   end
end
ndep = size(tt,1);
  

if ~exist([fnm '.nc'],'file')
   error(['S_FROM_DIRECT_TS05: Cannot find file ' fnm]);
end

rms1 = (nargout>=3);
rms2 = (nargout>=4);
   

X = getnc(fnm,'lon');
Y = getnc(fnm,'lat');
outofreg = find(lo<X(1) | lo>=X(end) | la<Y(1) | la>=Y(end));
if ~isempty(outofreg)
   disp([7 num2str(length(outofreg)) ' points outside region covered by t/s'])
   lo(outofreg) = [];
   la(outofreg) = [];
   tt(:,outofreg) = [];
   if ~isempty(doy)
      doy(outofreg) = [];
   end
end
ncast = length(lo);
ss = nan(size(tt));
if rms1
   rr1 = nan(size(tt));
end
if rms2
   rr2 = nan(size(tt));
end
if ncast==0
   return
end

% Begin to convert cast locations to indices into t-s climatology...
ix = interp1(X,1:length(X),lo);
iy = interp1(Y,1:length(Y),la);
ix0 = floor(min(ix));
ixe = floor(max(ix))+1;
iy0 = floor(min(iy));
iye = floor(max(iy))+1;
ix = 1+ix-ix0;
iy = 1+iy-iy0;
%X = X(ix0:ixe);
%Y = Y(iy0:iye);

ny = length(iy0:iye);

% Indices of neighbouring grid points (index = (nrows*(col-1)) + col)
i1 = round((ny*(floor(ix)-1))+floor(iy));
i2 = i1+1;
i3 = round((ny*floor(ix))+floor(iy));
i4 = i3+1;

% Calc horizontal and "vertical" interpolation weights
xr = ix-floor(ix);
yr = iy-floor(iy);
w = [(1-xr).*(1-yr); (1-xr).*yr; xr.*(1-yr); xr.*yr];

tlvs = getnc(fnm,'T_level');

tmp = find(tlvs<nanmin(tt(:)));
if isempty(tmp)
   itl(1) = 1;
else
   itl(1) = tmp(end);
end
tmp = find(tlvs>nanmax(tt(:)));
if isempty(tmp)
   itl(2) = length(tlvs);
else
   itl(2) = tmp(1);
end
tlvs = tlvs(itl(1):itl(2));

% Obtain S for T>32 by capping T at 32.   Justification:  any higher
% temperatures would normally be due to just diurnal heating, which wouldn't
% alter S.
ii = find(tt>=max(tlvs));
if ~isempty(ii)
   disp([num2str(length(ii)) ' T values capped for deriving S as they exceed' ...
	 ' range of S(T) climatology']);
   tt(ii) = max(tlvs)-.01;
end

timc = exp(doy*(-2i*pi/366));

MN = shiftdim(getnc(fnm,'mean',[itl(1) iy0 ix0],[itl(2) iye ixe]),1);  
if ann
   anc = shiftdim(getnc(fnm,'an_cos',[itl(1) iy0 ix0],[itl(2) iye ixe]),1);  
   AN = anc + ...
	1i.*shiftdim(getnc(fnm,'an_sin',[itl(1) iy0 ix0],[itl(2) iye ixe]),1);
   clear anc
end
if rms1
   R1 = shiftdim(getnc(fnm,'rmsr',[itl(1) iy0 ix0],[itl(2) iye ixe]),1);  
end
if rms2
   R2 = shiftdim(getnc(fnm,'rmsmr',[itl(1) iy0 ix0],[itl(2) iye ixe]),1);  
end

rr = find(isnan(MN));
MN(rr) = 0;
mn0 = squeeze(MN(:,:,1));
if ann
   rr = find(isnan(AN));
   AN(rr) = 0;
   an0 = squeeze(AN(:,:,1));
end   
if rms1
   R1(rr) = 0;
   r10 = squeeze(R1(:,:,1));
end   
if rms2
   R2(rr) = 0;
   r20 = squeeze(R2(:,:,1));
end   

for jj = 1:(length(tlvs)-1)
   mn = squeeze(MN(:,:,jj+1));   
   if ann
      an = squeeze(AN(:,:,jj+1));
   end
   if rms1
      r1 = squeeze(R1(:,:,jj+1));
   end
   if rms2
      r2 = squeeze(R2(:,:,jj+1));
   end
   ll = find(tt>=tlvs(jj) & tt<tlvs(jj+1));
   wz = (tt(ll)-tlvs(jj))./(tlvs(jj+1)-tlvs(jj));
   
   cst = ceil(ll/ndep);
   cst = cst(:)';
   
   if isempty(ll)
      ic = [];
   else
      % We cannot have a salinity of zero, so can use zero as the non-data flag.
      % If we left it as nan we would have to use the slower nansum below.
      indx = [i1(cst); i2(cst); i3(cst); i4(cst)];
      dd0 = mn0(indx);
      aa0 = ~~dd0;
      sumw0 = sum(aa0.*w(:,cst));   
   
      dd = mn(indx);
      aa = ~~dd;
      sumw = sum(aa.*w(:,cst));   
   
      % Require t/s values above and below each t, and that the data is not just 
      % at points almost the full grid interval away (ie that the good data 
      % interpolation weight is non-trivial)
      ic = find(sumw0>.05 & sumw>.05);
   end
   
   if ~isempty(ic)
      cic = cst(ic);
      
      aa0 = aa0(:,ic);
      dd0 = dd0(:,ic).*aa0;
      sumw0 = sumw0(ic);
      s1 = sum(dd0.*w(:,cic))./sumw0;

      if ann || rms1
	 indx = [i1(cic); i2(cic); i3(cic); i4(cic)];
      end
      
      if ann
	 dd0 = an0(indx).*aa0;
	 tmp = (sum(dd0.*w(:,cic))./sumw0);
	 s1 = s1 + real(tmp.*timc(cic));
	 an0 = an;
      end
      if rms1
	 dd0 = r10(indx).*aa0;
	 r1_1 = sum(dd0.*w(:,cic))./sumw0;
      end
      if rms2
	 dd0 = r20(indx).*aa0;
	 r2_1 = sum(dd0.*w(:,cic))./sumw0;
      end

      aa = aa(:,ic);
      dd = dd(:,ic).*aa;
      sumw = sumw(ic);
      s2 = sum(dd.*w(:,cic))./sumw;

      if ann
	 dd = an(indx).*aa;
	 tmp = sum(dd.*w(:,cic))./sumw;
	 s2 = s2 + real(tmp.*timc(cic));
      end
      if rms1
	 dd0 = r1(indx).*aa0;
	 r1_2 = sum(dd0.*w(:,cic))./sumw0;
	 rr1(ll(ic)) = r1_1'.*(1-wz(ic)) + r1_2'.*wz(ic); 
      end
      if rms2
	 dd0 = r2(indx).*aa0;
	 r2_2 = sum(dd0.*w(:,cic))./sumw0;
	 rr2(ll(ic)) = r2_1'.*(1-wz(ic)) + r2_2'.*wz(ic); 
      end

%      if size(s1,2)~=size(ic,2)
	 ss(ll(ic)) = s1'.*(1-wz(ic)) + s2'.*wz(ic);
%      else
%	 ss(ll(ic)) = s1.*(1-wz(ic)) + s2.*wz(ic);
%      end
   end
   
   mn0 = mn;
   if ann
      an0 = an;
   end
   if rms1
      r10 = r1;
   end
   if rms2
      r20 = r2;
   end
end
      
%--------------------------------------------------------------------------------
