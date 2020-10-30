% GET_ALL_CSL3  Get ocean cast data on CSIRO-standard-level v3 depths
%
%   NOTE:  No access yet to datasets 15  16  18  20  25-29  34
%
% INPUT:  
%  src    [w e s n] limits of region of required data
%      OR defining polygon [x1 y1; x2 y2; x3 y3; ... ] 
%  hvar   vector of one or more cast header info codes:
%         1)stnno   2)time   3)castflag  4)bot_depth   5)country   6)cruise
%         7)status [dset dependant]   9) filetype  11) profile depth  
%         12) bottom cast value  13) mixed layer depth
%            [3 11 & 12 apply to var(1) if multiple var]
%  var    vector of property codes: (see prop_name.m) 
%         1)t   2)s  3)02   4)Si  5)PO4   6)NO3   7)gamma  14)no4  15)nh3
%  dset   dataset codes - for info type:  dset_name
%  deps   Vector of indices of CSL v3 depth levels to extract
%         [use round(dep_csl(depths)) to convert from +ve metres to indices]
%         (Will interpolate to non CSL3 depths, but this is not recommended.)
%  scr    0 - disable pre-flagged bad-cast and bad-individual-data screening
%         1 - apply only originators' screening (eg NODC flags in WOD)
%         2 - [default] all originator & local screening
%  dups   0= don't remove any dups
%         1= only remove dups if "primary" datasets are also requested [default]
%         2= remove all "secondary" casts which have dups in "primary" datasets
%            (whether or not the primary datasets are being extracted now.)
%         dups{n} = [m o ...] =remove casts in n duped in datasets m,o,...
%         Use SHOW_DUP_PREFS to see from which datasets dup profiles are taken. 
%  strp   0= return all casts, even if no data
%         1= return only casts with some data in some profiles [default]
%         2= return only casts with some data in all requested properties
%  trng   [t1 t2] time range (days since 1900). [default - no restriction]
%
% GLOBAL variable SWITCHES:
%   Argo_Version   - deprecated, to select different versions of Argo archive
%   XBTtypesONLY    - 1=restrict QuOTA/IOTA to XBT probe types only
%   XBTFallRateCorr - 0=disable automatic fallrate correction of QuOTA XBTs
%
% OUTPUT:
%  lat,lon, then header vars [ncast 1] in order requested, then profile vars
%  [ncast ndep] in order {where ndep may be less than number of depths requested, 
%  if no good data that deep for that particular variable.}
%
% NOTE:  Empty depth columns will only occur if there is deeper good data.
%   If strp=0, all stations will be returned.
%   If strp=1, for each cast-row there will be at least some data in at least one
%   of the returned variables. 
%   If strp=2, there will be at least some data in all of the returned variables. 
%
% Q: Duplicates where more than one variable requested - eg request T & S,
%    one station has T only and the same station in another dataset has S only?   
% A: We have a different list of duplicates for each property, so in the
%    above case both stations would be returned (unless strp=2).
%
% USAGE: [lat,lon,v1,v2,..] = get_all_csl3(src,hvar,vars,dsets,deps,scr,dups,strp,trng);

% $Id: get_all_csl3.m,v 1.2 2006/01/04 05:55:53 dun216 Exp dun216 $
% Author: Jeff Dunn  CSIRO Marine Research Dec 1999
% Devolved from get_all_obs.m

function [lat,lon,varargout] = get_all_csl3(src,hvar,vars,dsets,deps,...
						scr,dups,strp,trng)

% Mods: 10/3/05 Extend dups system 
%       6/4/05  Add scr==2 code
%       23/8/06  Disabled auto-inclusion of dset 71 when ask for dset 7, var 6.
%       4/5/07  Global Argo_Version used to select between Argo dataset
%               versions
%       6/11/07 Wijffels et al 07 fall-rate correction to IOTA
%               Argo dmode available using hvar=7. 
%       12/2/08 Add MLD to hvars (so far only available in Argo)
%       26/5/08 Add IOTA Tasman (#37)
%       19/6/08 Add WOD05
%       26/4/11 Add WOD09
%       1/9/11  Add Pirata and Rama mooring arrays, update TAO
%       19/4/13 Cowley 2013 Fallrate and T-bias correction replaces Wijffels
%               correction for QuOTA XBTs.
%       9/5/13  Switches to filter out non-XBT QuOTA probes, and XBT fallrate
%               correction
%       4/5/17  Change french01_csl3 to french_csl3 (mon137)

global XBTtypesONLY
global XBTfallratecorWARN
persistent wod09_warned

if nargin<5 || isempty(deps)
   disp('  GET_ALL_CSL3  requires 5 or 6 input arguments')
   help get_all_csl3
   return
else
   deps = deps(:);
end

if nargin<6 || isempty(scr)
   scr = 2;
end   

if nargin<7 || isempty(dups)
   dups = 1;
end
if iscell(dups)
   % User has specified the dataset preference scheme for dup removal
   dupset = dups;
   if length(dupset)<max(dsets)
      dupset{max(dsets)} = [];
   end
   clear dups
   dups = 3;
elseif dups==1
   % Dup removal uses the standard preference scheme
   dupset = dsets;
elseif dups==2
   % Remove casts which duped in any other datasets, even if those datasets
   % are not being extracted	 
   dupset = [1:100];
end

if nargin<8 || isempty(strp)
   strp = 1;
end   

if nargin<9
   trng = [];
end   

if strp==0 & dups>0
   global Get_All_Csl3_told
   if ~Get_All_Csl3_told
      disp('You have set dups>0 & strp=0, so duplicate casts will have profiles')
      disp('emptied of data, but they will still be returned');
      Get_All_Csl3_told = 1;
   end
end

if max(deps>80) || any(diff(deps)==10) || any(deps==0)
   disp([7 'The "deps" vector given suggests you have specified depths in m,']);
   disp('rather than depth level indices! If nec., convert using dep_csl.');
end

mxdp = length(deps);
depm = csl_dep(deps,3);

ii = find(dsets==29);
if ~isempty(ii)
   pth = path_pc_or_nix('eez_data/boa_obslvl/wod01/SURF/surf_all.mat');
   disp('WOD01 SURF (dset 29) is not handled by GET_ALL_CSL. Access it')
   disp(['directly from ' pth]);
   dsets(ii) = [];
   if isempty(dsets)
      return
   end
end

ii = find(dsets<=6);
if ~isempty(ii)
   disp('There is no CSL3 version of WOD98 dsets (1:6)');
   dsets(ii) = [];
   if isempty(dsets)
      return
   end
end

if any(dsets>50 & dsets<=54) && isempty(wod09_warned)
   disp('      *     This warning once per session only      *')
   disp('WARNING: This version of WOD09 is a subset containing only profiles ')
   disp('with oxygen or nutrient data. Completely download again from NODC if')
   disp('wanting full set of T,S profiles.')
   disp('      *         *         *        *       *        *')
   disp(' ')
   wod09_warned = 1;
end
   
% 23/8/06 Disabled auto-inclusion of old CSIRO nitrate, because it appears
% to contain so much crap!
%if any(vars==6) & any(dsets==7) & ~any(dsets==71)
   % If want NO3 and CSIRO data, then also access the CSIRO Historical 
   % NO3 dataset (71). Added separately because this data was not in WOD98.
%   dsets = [dsets(:); 71];
%end
   
nhv = length(hvar);
ndv = length(vars);

lat = []; lon = [];

varargout{nhv+ndv} = [];

ndep = 0;

dpth = path_pc_or_nix('eez_data/boa_csl3/');
dupth = path_pc_or_nix('eez_data/boa_duplicates/');
infodir = path_pc_or_nix('eez_data/boa_qc_data/');
parnm = {'t','s','o2','si','po4','no3'};


for dset = dsets(:)'
   if any(dset==[21:28 41:48 51:54])
      [lo,la,stnno,hv,vv] = get_wod(dset,src,trng,deps,hvar,vars,scr,dpth);
   else
      [lo,la,stnno,hv,vv] = get_dset(dset,src,trng,deps,hvar,vars,scr,dpth);      
   end
   
   jj = 1:length(lo);
   
   if dups>0 
      % If dup checking, nan-fill all profiles in duplicate casts.
      % Messy if more than one var, because datasets have different
      % combinations of vars, so vars for the same cast can have
      % different dup status with respect to another dataset.
      % For each var we set dup profiles to nan. We then get rid of whole 
      % stations with only all-nan profiles.
      for jv = 1:ndv
	 if dups<3
	    dulst = dset_dup_pref(dupset,vars(jv),dset);
	 else
	    dulst = dupset{dset};
	 end
	 for odset = dulst
	    fnm = sprintf('%s%d_%d_%d_dups',dupth,dset,odset,vars(jv));
	    if exist([fnm '.mat'],'file')
	       load(fnm);
	       [tmp,ii] = intersect(stnno(jj),dupstn);
	       if ~isempty(ii)
		  vv{jv}(jj(ii),:) = nan;
	       end
	    end
	 end	 
      end   
   end

   if scr==2
      % Apply all local screening
      for jv = 1:ndv
	 if vars(jv)<=6 & ~isempty(vv{jv}) 
	    fnm = [infodir num2str(dset) filesep parnm{vars(jv)} '_scr'];
	    if exist([fnm '.mat'],'file')
	       load(fnm,'scrstn');
	       for jdep = 1:size(vv{jv},2)
		  II = find(ismember(stnno,scrstn{deps(jdep)}));
		  if ~isempty(II)
		     vv{jv}(II,jdep) = nan;
		  end
	       end
	    end
	 end
      end
   end
   
   % If removing casts with all (strp=1) or any (strp=2) empty profiles
   if strp>0
      if ndv==1
	 ii = find(all(isnan(vv{1}(jj,:)),2));
      else
	 bad = zeros([ndv length(jj)]);
	 for jv = 1:ndv
	    bad(jv,:) = all(isnan(vv{jv}(jj,:)),2)';
	 end
	 if strp==1
	    ii = find(all(bad));
	 else
	    ii = find(any(bad));
	 end
      end
      jj(ii) = [];
   end
   
   if ~isempty(jj)
      lat = [lat; la(jj)];
      lon = [lon; lo(jj)];   
      for ii = 1:nhv
	 varargout{ii} = [varargout{ii}; hv{ii}(jj)];
      end

      % For each variable...
      for ii = 1:ndv
	 % find depth of deepest data in any new cast
	 ldp = max(find(any(~isnan(vv{ii}(jj,:)),1)));
	 if isempty(ldp); ldp = 0; end	    
         [ncast,ndep] = size(varargout{nhv+ii});
	 
	 if ldp == ndep
	    % If existing output the same, just add to it
	    varargout{nhv+ii} = [varargout{nhv+ii}; vv{ii}(jj,1:ldp)];
	 elseif ldp < ndep
	    if size(vv{ii},2) >= ndep
	       varargout{nhv+ii} = [varargout{nhv+ii}; vv{ii}(jj,1:ndep)];
	    else
	       vpad = repmat(nan,[length(jj) ndep-ldp]);	       
	       varargout{nhv+ii} = [varargout{nhv+ii}; [vv{ii}(jj,1:ldp) vpad]];
	    end
	 else
	    % If new data deeper than existing output, pad out the output
	    vpad = repmat(nan,[ncast ldp-ndep]);
	    varargout{nhv+ii} = [[varargout{nhv+ii} vpad]; vv{ii}(jj,1:ldp)];
	 end
      end
   end
   
end                                 


return


%------------------------------------------------------------------------
function [lo,la,stn,hv,vv] = get_dset(dset,src,trng,deps,hvar,vars,scr,dpth);

global XBTtypesONLY
global XBTfallratecorWARN
global XBTFallRateCorr

if ~exist('XBTFallRateCorr','var') || isempty(XBTFallRateCorr)
   XBTFallRateCorr = 1;
end

la = [];
lo = [];
stn = [];
%stnno = [];
nhv = length(hvar);
ndv = length(vars);
hv{ndv} = [];
vv{ndv} = [];

% Header var names: 1-stn  2-time  3-cflag  4-botdep  5-co  6-cru  13-rejpc
%  9 is generated from 'dset'
% Header vars are all column vectors

inov = 1:ndv;
inoh = 1:nhv;
[iv,ih] = dset_vars(dset,2,vars,hvar);
if isempty(iv)
   return
end
inov(iv) = [];
inoh(ih) = [];


vnm = {'t','s','o2','si','po4','no3','nutdens','','','','','','','no4','nh3'};

cflag = [];
for ii = 1:ndv
   eval([vnm{vars(ii)} ' = [];']);
   eval([vnm{vars(ii)} '_castflag = [];']);
end

lon = []; lat = [];

if any(dset==[8 17 37 38])
   if isempty(XBTfallratecorWARN)
       XBTfallratecorWARN = 1;
       disp(['NOTE: The QuOTA XBT data via GET_ALL_CSL3 is now Cowley 2013' ...
	     ' fallrate corrected. WOD XBT datasets are not presently corrected' ...
	     '- ask Jeff if this is required']);       
       disp('This notification provided once per session.');
   end
end

wantcflg = ((scr>0  && vars(ii)<=6) || any(hvar==3));

switch dset
   
  case 7
    % CSIRO
    fnm = [dpth 'csiro' filesep 'csiro_csl'];
    lat = getnc(fnm,'lat');
    lon = getnc(fnm,'lon');
    time = getnc(fnm,'time');
    stnno = getnc(fnm,'stnno');
    if any(hvar==4)
       botdep = getnc(fnm,'botdepth');
    end
    if any(hvar==5)
       co = repmat(12345,size(lat));
    end
    if any(hvar==6)
       cru = floor(stnno/1000);
    end
    %if any(hvar==11 | hvar==12)
    %   bnm = [vnm{vars(1)} 'bot'];
    %   eval([bnm ' = getnc(fnm,bnm);']);
    %end
    for ii = iv
       nm = vnm{vars(ii)};
       eval([nm ' = getnc(fnm,nm);']);
       if wantcflg
	  flnm = [nm '_castflag'];
	  eval([flnm ' = getnc(fnm,flnm);']);
       end    
    end
    
  case 8
    % IOTA East
    
    % Calc and restrict stdep as only 55 levels in IOTA file.
    deps = deps(find(deps<=55));

    load([dpth 'csiro_therm_archive' filesep 'iota_east_csl3']);

    botdep = bdep;
    t = tz;    
    %if ~isempty(XBTtypesONLY) && XBTtypesONLY
       % Known XBT types, or XBT of unknown type but Haniwa-corrected  
       % Could also use dtyp==22594 & probe_type==0 : ie XBTs of unknown type
     %  disp('Reducing IOTA to just XBTs');
     %  jj = ((probe_type>0 & probe_type<700) | probe_type==992);
     %  t = t(jj,:);
     %  botdep = botdep(jj);
     %  lat = lat(jj);
     %  lon = lon(jj);
     %  time = time(jj);
     %  stnno = stnno(jj);
     %  deep = deep(jj);
     %  TSK = TSK(jj);
     %  Tnum = Tnum(jj);
     %  t_castflag = t_castflag(jj);
    %end
    clear tz bdep dtyp indx

    
  case 9
    % NIWA
    load([dpth 'niwa_csl3']);
    cru = floor(stnno/10000);
    if any(hvar==5)
       co = repmat(13873,[length(lon) 1]);
    end
    
  case 10 
    % French Indian 2001 CD
    load([dpth 'french_csl3'])
    
    % These cruise codes are crap - almost worth not using them!
    cru = cru_num;    
    
  case 11
    % ARGO floats  (only available to depth 66 (2000m))
    deps = deps(find(deps<=66));
    global Argo_Version Argo_V_Warn1 Argo_V_Warn2
    if isempty(Argo_Version)
       if isempty(Argo_V_Warn1) & exist([dpth 'Argo_Update_Active.mat'],'file')
	  load([dpth 'Argo_Update_Active'],'note');
	  disp('****')
	  disp(note)
	  Argo_V_Warn1 = 1;
       end
       load([dpth 'argo_csl3'])
    else
       vdesc = {'7/9/07 Pressure-corrected','3/5/10','latest'};
       Avnm = {'argoPC_csl3','argo_csl3_03May2010','argo_csl3'};
       av = Argo_Version;
       if isempty(Argo_V_Warn2)
	  if ~any([1 2 3]==av)
	     disp('WARNING - GET_ALL_CSL3: not a valid value for Argo_Version')
	     disp('          Setting Argo_Version=3 for now')
	  else
	     disp(['WARNING - using ' vdesc{av} ' version of Argo'])
	     if scr>1
		disp('and any per-value screening may mismatch with this dataset');
	     end
	  end
	  Argo_V_Warn2 = 1;
       end

       load([dpth Avnm{av}]);
    end
    if isa(deps,'single')
       deps = double(deps);
    end
    
  case {12,121}
    % TAO moorings  (ignores tiny amount of data available to depth 66 (2000m))
    deps = deps(find(deps<=55));    
    if isempty(deps)
       return
    end
    
    % tao_csl3 has 333870 daily profiles. For mapping purposes, efficient to
    % instead use tao_month_csl3 (dset 121), the seasonal averages! 
    % 'cru' is our consecutive mooring number.
    if dset==12
       load([dpth 'tao_csl3']);
    else
       load([dpth 'tao_month_csl3']);
    end
  
  case 13
    % Antarctic CRC
    fnm = [dpth 'crc_csl3'];
    lat = getnc(fnm,'lat');
    lon = getnc(fnm,'lon');
    time = getnc(fnm,'time');
    stnno = getnc(fnm,'stnno');
    if any(hvar==4)
       botdep = getnc(fnm,'botdepth');
    end
    if any(hvar==5)
       co = repmat(12345,[length(lat) 1]);
    end
    if any(hvar==6)
       % cru = getnc(fnm,'cru_idx');  % Replaced 23/5/08 because not useful
       cru = floor(stnno/1000);
    end
    for ii = iv
       nm = vnm{vars(ii)};
       eval([nm ' = getnc(fnm,nm);']);
       if wantcflg
	  flnm = [nm '_castflag'];
	  eval([flnm ' = getnc(fnm,flnm);']);
       end
    end    
       
       
  case 14
    % Willis global QC-ed thermal
    % This is a 750MB dataset, so is split into chunks. get_willis collates it.
    % Also because it is huge, we select only the required depths. To make
    % this work, need to adjust "deps" so that it corresponds to the index of
    % the extracted data. eg if want depths 50 60 75 120, ie Willis depths 
    % [5 6 7.5 12], then extract 5:12 and adjust "deps" to [1 2 3.5 8].
    
    % Convert CSLv3 indices to Willis depth indices (non-integer is an (n+5)m
    % CSL depth required).
    wdep = 0:10:750;
    deps = interp1(wdep,1:length(wdep),csl_dep(deps,3));
    deps = deps(find(deps>=1 & deps<=length(wdep)));
    
    if isempty(deps)
       return
    end
       
    % Reduce to minimum required block of depths, and adjust index to start
    % at 1 for the first of these depths.
    ii = min(floor(deps)):max(ceil(deps));
    deps = deps+1-ii(1);
    wdep = wdep(ii);

    % Note that get_willis returns non-TAO Willis by default. 
    tao14 = 0;
        
    [lon,lat,time,t,botdep,a1,a2,a3,a4,stnno] = get_willis(src,wdep,tao14);       
    clear a?
        
  case 15
    % #### French Pacific SSS
    disp('No access to Pacific SSS')

  case 16
    % #### Far Seas
    disp('No access to Far Seas yet, but is only CSLv1 anyway')

  case 17
    % IOTA West XBT
    % Calc and restrict stdep as only 55 levels in IOTA file.
    deps = deps(find(deps<=55));

    load([dpth 'csiro_therm_archive' filesep 'iota_west_csl3']);
    botdep = bdep;
    t = double(tz);

    stnno = double(stnno);
    
    clear tz bdep dtyp indx
        
  case 18
    % #### MEDS GTS XBT etc
    disp('No access to MEDS GTS QCed Thermal yet!')

  case 19
    % WOCE WHP
    load([dpth 'woce_csl3']);
    if any(hvar==6)
       % Any WOCE station may be associated with up to 3 cruise codes
       % (which were used to make up the original netcdf station file name.)
       % The 352 unique codes are in variable 'cruz'. 'cri1' points to the 
       % first cruise code for each station, 'cri2' to the second (if
       % required), 'cri3' to the third, if required. We can't pass that 
       % complication through this general function, so just use cri1.
       cru = cri1;
    end
    
  case 20
    % #### WOCE UOT
    disp('No access to WOCE UOT yet!')
  
  case 31
    % Other small dsets, esp AIMS Torres CTD 
    % Restrict deps as only 11 levels in Torres file.
    deps = deps(find(deps<=11));

    load([dpth 'other_csl3']);   

    if any(hvar==5)
       co = repmat(12345,[length(lon) 1]);
    end
    
  case 32
    % Bernadette Sloyan SE Pacific
    load([dpth 'sloyan_SEPac_csl3']);   

  case 34
    % #### High density XBT
    disp('No access to dset 34 [high density XBT] yet!')
  
  case 35
    % AWI (Alfred Wegener Institute) CTD
    load([dpth 'awi_csl3']);    
  
  case 37
    % IOTA Tasman    
    % Calc and restrict stdep as only 55 levels in IOTA file.
    deps = deps(find(deps<=55));
    
    load([dpth 'csiro_therm_archive' filesep 'iota_tasman_csl3']);
    botdep = bdep;
    t = double(tz);
    
    clear tz bdep dtyp indx
    
  case 38
    % Combined IOTA - all probe types
    % Calc and restrict stdep as only 55 levels in IOTA file.
    deps = deps(find(deps<=55));

    load([dpth 'csiro_therm_archive' filesep 'iota_all_csl3']);
    botdep = bdep;
    t = double(tz);

    stnno = double(stnno);
    
    clear tz bdep dtyp indx probe_type maxdep
        
  case {61,161}
    % Pirata moorings 
    deps = deps(find(deps<=45));    
    if isempty(deps)
       return
    end
    
    % For mapping purposes, efficient to instead use monthly avs (dset 161)
    if dset==61
       load([dpth 'pirata_csl3']);
    else
       load([dpth 'pirata_month_csl3']);
    end
  
  case {62,162}
    % RAMA moorings  
    deps = deps(find(deps<=50));    
    if isempty(deps)
       return
    end
    
    % For mapping purposes, efficient to instead use monthly avs (dset 162)
    if dset==62
       load([dpth 'rama_csl3']);
    else
       load([dpth 'rama_month_csl3']);
    end
  
  case 70
    disp('Dataset 70 (CSIRO 2db CTD) is an obs-level-only dataset.')
  
  case 71
    % CSIRO historical nitrate
    fnm = [dpth 'csiro' filesep 'csiro_no3'];
    lat = getnc(fnm,'lat');
    lon = getnc(fnm,'lon');
    time = getnc(fnm,'time');
    stnno = getnc(fnm,'stnno');
    if any(hvar==4)
       botdep = getnc(fnm,'botdepth');
    end
    if any(hvar==6)
       cru = floor(stnno/1000);
    end
    for ii = iv
       nm = vnm{vars(ii)};
       eval([nm ' = getnc(fnm,nm);']);
       if wantcflg
	  flnm = [nm '_castflag'];
	  eval([flnm ' = getnc(fnm,flnm);']);
       end
    end    

  case 72
    disp('Dataset 72 (CSIRO Hydro) is an obs-level-only dataset.')
  
  case 81
    % ENACT3 WOD05 CTDs 1997-2004  Temporary
    load([dpth 'wod05_en3_ctd_csl3'])
    stnno = zeros(size(lon));
    cru = zeros(size(lon));
    
  otherwise
    disp(['Do not yet have access to dataset ' num2str(dset)]);
    
end

if exist('stnno','var'); stnno = double(stnno); end
if exist('botdep','var'); botdep = double(botdep); end
if exist('cru','var'); cru = double(cru); end
    

if isempty(src)
   inreg = 1:length(lon);
elseif min(size(src))==1
   inreg = find(lon>=src(1) & lon<=src(2) & lat>=src(3) & lat<=src(4));
else
   inreg = find(inpolygon(lon,lat,src(:,1),src(:,2)));
end
    
if ~isempty(trng) & ~isempty(inreg)
   kk = find(time(inreg)<trng(1) | time(inreg)>=trng(2));
   inreg(kk) = [];
end

if ~isempty(inreg)
   % When loading all the vector variables, we ensure they are columns.
   nin = length(inreg);
       
   if (dset==8 || dset==17 || dset==37 || dset==38)
      if ~isempty(XBTtypesONLY) && XBTtypesONLY
	 % Known XBT types, or XBT of unknown type but Haniwa-corrected  
	 % Could also use dtyp==22594 & probe_type==0 : ie XBTs of unknown type
	 disp('Reducing IOTA to just XBTs');
	 jj = ((probe_type(inreg)>0 & probe_type(inreg)<700) | probe_type(inreg)==992);
	 inreg = inreg(jj);
      end
      
      if XBTFallRateCorr 
	 % Wijffels 2008 correction
	 %t(inreg,:) = xbt_fallrate_corr(t(inreg,:),z,deep(inreg),time(inreg),lo,la);
	 % Cowley 2013 correction
	 t(inreg,:) = xbt_fallrate_corr_13(t(inreg,:),z,deep(inreg),time(inreg),TSK(inreg),Tnum(inreg));
      end
   end

   la = lat(inreg);      la = la(:);
   lo = lon(inreg);      lo = lo(:);
   stn = stnno(inreg);   stn = stn(:);

   for ii = ih
      switch hvar(ii)
	case 1
	  hv{ii} = stn;
	case 2
	  hv{ii} = time(inreg);
	case 3
	  if isempty(cflag) & ~isempty(iv)
	     cfnm = [vnm{vars(iv(1))} '_castflag'];
	     if ~isempty(eval(cfnm))
		eval(['cflag = ' cfnm ';']);
	     end
	  end
	  if ~isempty(cflag)
	     hv{ii} = cflag(inreg);
	  end
	case 4
	  hv{ii} = botdep(inreg);
	case 5
	  hv{ii} = co(inreg);
	case 6
	  hv{ii} = cru(inreg);
	case 7
	  if dset==11
	     hv{ii} = dmode(inreg);
	  else
	     hv{ii} = repmat(nan,[nin 1]);
	  end
	case 9
	  hv{ii} = repmat(dset,[nin 1]);
	case 11
	  hv{ii} = eval([vnm{vars(1)} 'bot(inreg,1)']);
	case 12
	  hv{ii} = eval([vnm{vars(1)} 'bot(inreg,2)']);
	case 13
	  hv{ii} = double(rejpc(inreg));
	otherwise
	  hv{ii} = repmat(nan,[nin 1]);		 
      end
      hv{ii} = hv{ii}(:);
   end
   
   for ii = inoh
      hv{ii} = repmat(nan,[nin 1]);
   end	

   ndp = length(deps);
   
   for ii = iv
      if isempty(eval([vnm{vars(ii)} '_castflag'])) & ~isempty(cflag)
	 % One castflag for all variables, so use an approp named copy
	 eval([vnm{vars(ii)} '_castflag = cflag;']);
      end
      switch vars(ii)
	case 1
	  vv{ii} = [vv{ii}; scrload(t,inreg,deps,t_castflag,scr)];
	case 2
	  vv{ii} = [vv{ii}; scrload(s,inreg,deps,s_castflag,scr)];
	case 3
	  vv{ii} = [vv{ii}; scrload(o2,inreg,deps,o2_castflag,scr)];
	case 4
	  vv{ii} = [vv{ii}; scrload(si,inreg,deps,si_castflag,scr)];
	case 5
	  vv{ii} = [vv{ii}; scrload(po4,inreg,deps,po4_castflag,scr)];
	case 6
	  vv{ii} = [vv{ii}; scrload(no3,inreg,deps,no3_castflag,scr)];
	case 7
	  vv{ii} = nutdens(inreg,deps);
	case 14
	  vv{ii} = [vv{ii}; scrload(no4,inreg,deps,no4_castflag,scr)];
	case 15
	  vv{ii} = [vv{ii}; scrload(nh3,inreg,deps,nh3_castflag,scr)];
	otherwise
	  vv{ii} = repmat(nan,[nin ndp]);		 
      end
   end

   for ii = inov
      vv{ii} = nan(nin,ndp);
   end
end


return

%------------------------------------------------------------------------
function vo = scrload(vv,jj,deps,cflg,scr)

% Are some of the deps BETWEEN rather than ON CSL3 levels?
intp = any(round(deps)~=deps);

if intp
   % If the data depth is close to a CSL depth, we say that is close enough
   % (these are indices, so .1 refers to portion of depth *interval*)
   if max(abs(round(deps)-deps))<.1
      deps = round(deps);
      intp = 0;
   else
      % Use 1q because it doesn't stop with a warning if nan's in data.
      vo = interp1q([1:size(vv,2)]',vv(jj,:)',deps)';
   end
end
if ~intp
   vo = vv(jj,deps);
end

if scr & ~isempty(cflg)
   kk = find(cflg(jj)>0);
   vo(kk,:) = nan;
end

return

%------------------------------------------------------------------------
% WOD is different because stored in WMO files.

function [lo,la,stn,hv,vv] = get_wod(dset,src,trng,deps,hvar,vars,scr,dpth);


la = [];
lo = [];
stn = [];
nhv = length(hvar);
ndv = length(vars);
hv{nhv} = [];
vv{ndv} = [];

if any(dset==[25 26 27 28 29 45 47 48 49 55 56 57 58])
   disp(['No access to dataset ' num2str(dset) ' yet!'])
   return
end
  
vnm = {'t','s','o2','si','po4','no3','nutdens'};

wod05 = (dset>40 & dset<=50);
wod09 = (dset>50 & dset<=59);

if wod09
   pth = [dpth 'wod09' filesep];
   prefx = {'CTD','OSD','PFL','UOR','DRB','GLD','MRB','XBT'};
   ids = dset-50;
elseif wod05
   pth = [dpth 'wod05' filesep];
   prefx = {'CTD','OSD','PFL','UOR','DRB','GLD','MRB','XBT'};
   ids = dset-40;
else
   pth = [dpth 'wod01' filesep];
   prefx = {'CTD','OSD','PFL','UOR','DRB','MBT','MRB','XBT'};
   ids = dset-20;
end

if isempty(src)
   if wod09
      src = [0 360 -80 90];
   else
      src = [0 360 -80 90];
   end      
end    
if min(size(src))==1
   wmosq = getwmo(src);
else
   wmosq = getwmo(...
       [min(src(:,1)) max(src(:,1)) min(src(:,2)) max(src(:,2))]);
end


for wmo = wmosq   
   fnm = [pth prefx{ids} filesep prefx{ids} num2str(wmo) 'csl.mat'];

   if exist(fnm,'file')
      F = load(fnm);
      if min(size(src))==1
	 jj = find(F.lon>=src(1) & F.lon<=src(2) & F.lat>=src(3) & F.lat<=src(4));
      else
	 jj = find(inpolygon(F.lon,F.lat,src(:,1),src(:,2)))';
      end

      if ~isempty(trng) & ~isempty(jj)
	 kk = find(F.time(jj)<trng(1) | F.time(jj)>=trng(2));
	 jj(kk) = [];
      end

      if ~isempty(jj)
	 nj = length(jj);

	 la = [la; F.lat(jj)];
	 lo = [lo; F.lon(jj)];
	 stn = [stn; F.ocl(jj)];
	 
	 for ii = 1:nhv
	    switch hvar(ii)
	      case 1
		hv{ii} = [hv{ii}; F.ocl(jj)];
	      case 2
		hv{ii} = [hv{ii}; F.time(jj)];
	      case 3
		if ~isempty(vars)
		   gotflg = 0;
		   cfnm = [vnm{vars(1)} '_castflag'];
		   if isfield(F,cfnm) 
		      cflag = getfield(F,cfnm);
		      if ~isempty(cflag)
			 hv{ii} = [hv{ii}; cflag(jj)];
			 gotflg = 1;
		      end
		   end
		   if ~gotflg
		      hv{ii} = [hv{ii}; zeros(length(jj),1)];
		   end
		end
	      case 4
		hv{ii} = [hv{ii}; F.botdep(jj)];
	      case 5
		hv{ii} = [hv{ii}; F.cc(jj)];
	      case 6
		hv{ii} = [hv{ii}; F.cru(jj)];
	      case 9
		hv{ii} = [hv{ii}; repmat(dset,[nj 1])];
	      case {11,12}
		gotb = 0;
		vfnm = [vnm{vars(1)} 'bot'];
		if isfield(F,vfnm) 
		   botv = getfield(F,vfnm);
		   if ~isempty(botv)
		      gotb = 1;
		      if hvar(ii)==1
			 hv{ii} = [hv{ii}; botv(jj,1)];
		      else
			 hv{ii} = [hv{ii}; botv(jj,2)];
		      end
		   end
		end
		if ~gotb
		   hv{ii} = [hv{ii}; repmat(nan,[length(jj) 1])];
		end
	      otherwise
		hv{ii} = [hv{ii}; repmat(nan,[nj 1])];		 
	    end
	 end	

	 for ii = 1:ndv
	    if vars(ii)<=7 & isfield(F,vnm{vars(ii)})
	       switch vars(ii)
		 case 1
		   vv{ii} = [vv{ii}; scrload(F.t,jj,deps,F.t_castflag,scr)];
		 case 2
		   vv{ii} = [vv{ii}; scrload(F.s,jj,deps,F.s_castflag,scr)];
		 case 3
		   vv{ii} = [vv{ii}; scrload(F.o2,jj,deps,F.o2_castflag,scr)];
		 case 4
		   vv{ii} = [vv{ii}; scrload(F.si,jj,deps,F.si_castflag,scr)];
		 case 5
		   vv{ii} = [vv{ii}; scrload(F.po4,jj,deps,F.po4_castflag,scr)];
		 case 6
		   vv{ii} = [vv{ii}; scrload(F.no3,jj,deps,F.no3_castflag,scr)];
		 case 7
		   vv{ii} = [vv{ii}; F.nutdens(jj,deps)];
		 otherwise
		   vv{ii} = [vv{ii}; repmat(nan,[nj length(deps)])];
	       end
	    else  
	       vv{ii} = [vv{ii}; repmat(nan,[nj length(deps)])];
	    end
	 end
      end    % endif ~isempty(jj)
   
      clear F
   end    % endif exist(fnm,'file')
end     % endfor wmo = wmosq 
	 

%-----------------------------------------------------------------------------
