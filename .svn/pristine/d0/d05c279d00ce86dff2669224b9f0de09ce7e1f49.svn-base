% DEPS_TO_STDDEP  Interpolate values on arbitrary depths to standard depths
%   (Designed for bottle data, but fine for 2db CTD data because it handles
%    any gaps safely, although the rr interpolation is a waste of time for
%    CTDs as it would very rarely be applicable or beneficial - so for CTDs
%    use maxdis=-1 )
%
%    NOTE: least squares linear extrapolation used over small distance if no 
%      interp done (esp beyond range of obs). Errors caused by using this 
%      can be reduced using the "grad_lim" and "botlim" options, and with 
%      noisy high-res data using the "noisy" option (see below).
%        
% INPUT: 
%  odep  - *vector* of depths of the data (in m)
%  obs   - corresponding data values, with NaNs indicating any gaps.
%    OPTIONAL:
%  maxdis - see doco in "help rr_int"   [default 1 - use default spacing]
%           A value of -1 mean this code will NOT use the rr_int system.
%  sdep   -  depth levels to use. Defaults to Standard Depths
%           Depths must be +ve and increasing.
%  noisy -  (suitable only for high-res potentially noisy data)
%           1= extrapolate by least squares linear fit  [default]
%           2= pre-smooth by 2m bin mean filter and least squares extrap
%           3= locally-weighted-mean filter to 2m and least squares extrap
%            [2 & 3 are good where noise > dv/dz, and bad where dv/dz >> noise.]
%  botlmt - 1=do NOT extrapolate much below the deepest
%           observations [default 0]
%  grad_lim - if given, will not extrapolate if it would exceed the specified  
%           gradient (dV/dz) eg PSU/m, or total offset exceeds
%           10*grad_lim. This should of course be property-specific, but no 
%           facility yet to be depth or location specific.  
%           eg could use : T .25 degC/m, => total offset <2.5 degC
%           suggest  T .25,   S .1 to .05,  O2 .1,   PO4 .025,  SiO&NO3 .25
%           ALSO note that extrapolation can create imposisble values, eg
%           -ve nutrients. Cannot automatically test for this (unless we 
%           pass in limit values) so need to check afterwards.
% 
% OUTPUT: sdat - interpolated values on standard depths levels
%
%    Also, global variable activity counters are updated
%
% MODS:  
% 16/6/05 remove old linear extrap code and now always use least-squares 
%         linear. Also only direct-subs if can't extrap, and restrict extrap
%         if obs span is less than extrap distance.
%
% $Id: deps_to_stddep.m,v 1.11 2005/06/16 02:18:35 dun216 Exp dun216 $
% Jeff Dunn  5/5/97  Copyright CSIRO Marine Research
%
% USAGE:  sdat = deps_to_stddep(odep,obs)
%    OR:  sdat = deps_to_stddep(odep,obs,[],[0 10 50 1000 1010])
%    OR:  sdat = deps_to_stddep(odep,obs,maxdis)
%           after first call to get threshold matrix from rr_int:
%             [dum,maxdis] = rr_int([]);

function [sdat] = deps_to_stddep(odep,obs,maxdis,sdep,noisy,botlmt,grad_lim)

% Previous MODS Notes:
% 11/6/98  remove depth inversions
%
% 22/7/02 Mods and tuning (esp the extrapolation code) to suit the salinity
%   data of the French CD of 2001, where many crappy pre-smoothed 10db casts. 
%   (changed near_lim etc from step wise to linear functions, and extrapolate 
%   instead of most cases of direct substitution.)
%
% 14/12/04 Introduce options noisy=2 & 3, to pre-filter hi-res casts; option 
%   grad_lmt; and refine botlmt.
%
% 11/5/2011 Place extra restriction on selection of linear interpolation
%   cases to prevent interp across too large a gap.


% Global cummulative activity counters:
% rr_int_cnt  - rr_int interpolated values
% lin_int_cnt - linear interpolated
% dir_sub_cnt - direct substitutoin
% l_extrp_cnt - linear extrapolation (from 1 or 2 pairs of points)
% r_extrp_cnt - regression-fit extrapolation of up to 6 potentially noisy
%               points (combined with linear extrap of nearest pair.)

global rr_int_cnt lin_int_cnt dir_sub_cnt l_extrp_cnt r_extrp_cnt;

if nargin<7
   grad_lim = [];
else
   offlim = 10*grad_lim;
end

if nargin<6 | isempty(botlmt)
   botlmt = 0;
end

if nargin<5 | isempty(noisy)
   noisy = 0;
end

if nargin<4 | isempty(sdep)
   nlvl = dep_csl(5000,3);
   sdep = csl_dep(1:nlvl,3)';
end

odep = odep(:);
obs = obs(:);
sdep = sdep(:);
nlvl = length(sdep);

xfn = [0 300 1200 8000];
yfn = [7 15 75 150];
near_lim = interp1(xfn,yfn,sdep);
far_lim = 2*near_lim;
dir_lim = near_lim/5;
% => dir_lim direct substitution limit = 1.4 @ 0, 15 @ 1200, 30 @ 8000m

sdat = repmat(NaN,nlvl,1);

% Remove NaNs, test depths range (as DPG files are sometimes crap), and remove
% depth inversions.

jj = find(isnan(obs) | isnan(odep) | odep<0 | odep>8000);
if ~isempty(jj)
  obs(jj) = [];
  odep(jj) = [];
end
if ~isempty(obs)
  jj = find((odep(2:end)-odep(1:end-1))<=0);
  if ~isempty(jj)
    obs(jj+1) = [];
    odep(jj+1) = [];
  end
end
ndeps = length(obs);

if ndeps == 0
   % RETURN if no data
  return;
end


% Pre-smoothing
if noisy>1
   vonly = 1;
   if noisy==1
      lnscl = 1.05;   % approx bin means at 2m spacing
   else
      lnscl = 3;
   end
   xi = 0:2:max(odep);
   [obs,odep] = filt_y_to_x(odep,obs,xi,noisy,vonly,lnscl);
   if isempty(obs)
      return   
   end
   ndeps = length(obs);
end

if nargin<3 | isempty(maxdis); maxdis = 1; end
   
if ndeps < 4 | maxdis == -1
  sidx = (1:nlvl)';
else  
   % RR INTERPOLATION
   sdat = rr_int(obs,odep,sdep,1,maxdis);
   sidx = find(isnan(sdat));
   rr_int_cnt = rr_int_cnt + nlvl - length(sidx);
end


% LINEAR INTERPOLATION
% Find which depths left to fill by linear interpolation. First, get index
% to sdep's between the odep's, then find between which pairs of odep's each
% remaining sdep falls, then calc distance from the two odep's, and assess 
% whether either near enough to one odep or the gap between the two odep's is 
% small enough to interpolate a value at that sdep.

if ~isempty(sidx)  & ndeps >= 2
  idx = sidx(find(sdep(sidx)>odep(1) & sdep(sidx)<odep(ndeps)));

  if ~isempty(idx)
    oidx = interp1(odep,1:ndeps,sdep(idx));
    dists = [sdep(idx)-odep(floor(oidx)) odep(ceil(oidx))-sdep(idx)];  
    near = min(dists')';
    far = max(dists')';

    %# interp = idx(find(near<near_lim(idx) | far<far_lim(idx)));
    % Refinement to cut out more extreme linear interp  11/5/2011
    interp = idx(find(far<far_lim(idx) | (near<near_lim(idx) & far<2*far_lim(idx))));

    if ~isempty(interp)
      sdat(interp) = interp1(odep,obs,sdep(interp));
      sidx = find(isnan(sdat));
      
      lin_int_cnt = lin_int_cnt + length(interp);
    end
  end
end


% DIRECT SUBSTITUTION & EXTRAPOLATION
% If values still unfilled, see if almost or exactly the same depth as an obs
% so can just take that value (ie direct substitution), or otherwise if near
% enough to at least 2 obs so can safely linearly extrapolate. 

% Should have a gradient limit on extrapolation? This would need to be
% property-specific (and maybe depth-specific), so is a bit tricky!

if ~isempty(sidx)
   % idx is nearest obs to the remaining sdeps. Any sdep beyond the range of 
   % odep would crash interp1 - prevent by adding extreme values to odep and
   % indexing these 0 and ndep+1 (so that odep are still 1 to ndep).
   
   idx = round(interp1([-99999; odep; 99999],0:ndeps+1,sdep(sidx)));
 
   kk = find(abs(odep(idx)-sdep(sidx)) < near_lim(sidx));
   for jj = kk(:)'
      sdj = sdep(sidx(jj));
      odj = odep(idx(jj));
      x = sdj-odj;      
      new = nan;
      
      if abs(x) > 1.5
	 % Only extrapolate if last obs is more than 1.5m from target depth
	 jll = find(abs(odep-sdj) < far_lim(sidx(jj))); 
	 if x > 0
	    jll = flipud(jll);
	 end
	 
	 if length(jll)<2 | max(abs(odep(jll)-odj)) < abs(x)
	    jll = [];
	 elseif any(abs(diff(odep(jll))) < 1.5)
	    ll = jll(1);
	    for mm = jll(2:end)'
	       if abs(odep(ll(end))-odep(mm))>1.5  & ...
	       (length(ll) < 4 | abs(odj-odep(mm)) < abs(x))
		  ll = [ll mm];
	       end
	    end
	    jll = ll;
	 end
            
	 if length(jll) >= 2
	    if abs(max(obs(jll))-min(obs(jll)))<.005
	       % If obs basically constant, just use nearest (as it will
               % crash the line-fit below anyway)	       
	       new = obs(jll(1));
	    else
	       % line-fit algorithm corrected on 18/6/08
	       xog = min(odep(jll));
	       cc = ([ones([length(jll) 1]) odep(jll)-xog]) \ obs(jll);
	       new = cc(1) + (sdj-xog)*cc(2);
	    end
	    r_extrp_cnt = r_extrp_cnt + 1;
	    
	    if ~isempty(grad_lim)
	       ofset = abs(obs(idx(jj))-new);
	       if ofset>abs(x*grad_lim) | ofset>offlim 
		  new = nan;
	       end
	    end
	    sdat(sidx(jj)) = new; 
	 end
      end
      
      if isnan(new) & abs(x)<dir_lim(sidx(jj))
	 % Couldn't extrapolate, but within direct subs range, so do it!
	 sdat(sidx(jj)) = obs(idx(jj)); 
	 dir_sub_cnt = dir_sub_cnt + 1;
      end
   end
   sidx = find(isnan(sdat));
end

if botlmt
   bval = sdep>(odep(ndeps)+10);
   if any(bval)
      sdat(bval) = nan;
   end

   % Pre 2011:  reject at most just the bottom-most value
   %bval = max(find(~isnan(sdat)));
   %if ~isempty(bval) & (sdep(bval)-odep(ndeps))>10
   %   sdat(bval) = nan;
   %end
end
   
   
% ------------ End of deps_to_stddep -----------------
