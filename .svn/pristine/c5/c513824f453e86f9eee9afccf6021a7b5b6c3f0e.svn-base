% GET_ALL_OBS  Get BOA ocean cast data on Observed levels
%
%   NOTE: This superceded and is incompatible with previous version, which 
%        is still available but now called get_all_obs_v1.
%
%   NOTE:  No access yet to datasets 8  15  16  17  18  20  25-29  32 34
%
% INPUT:  
%  src    [w e s n] limits of region of required data
%      OR defining polygon [x1 y1; x2 y2; x3 y3; ... ] 
%  hvar   vector of *ONE or MORE* cast header info codes:
%         1)stnno   2)time   3)castflag  4)bot_depth   5)country   6)cruise
%         7)DMode 9) filetype 13)rejpc  [NB castflag is for first specified var]
%  var    vector of property codes:  
%         1)t   2)s   3)02   4)Si   5)PO4   6)NO3   7)gamma 9)pH  14)no4  15)nh3
%  dset   dataset codes - for info type:  dset_name
%         WARNING: some dsets, especially #14, can use enormous memory. #14
%         may be better accessed directly using get_willis.m
%  dlim   [Optional] upper and lower depth limits
%  dmid   [Optional] target depth (return nearest nval values to this depth)
%  nval   [Optional] max number of values to return (in dlim range, closest
%          dmid if dmid given, otherwise limited trim off at bottom of
%          profile. [default 1 if dmid, otherwise 25000] 
%  scr   -1 - disable pre-flagged bad-cast and bad-individual-data screening
%         0 - only     "      "        "        "     "      "       "
%         1 - apply castflag screening
%         2 - [default] apply (originator and CSIRO) per-value screening as well
%         3 - screen as for 2, but leave all rejects in casts [used when doing
%             further screening, so that value indices match original casts.]
%  trng   [mint maxt] optional time limits (days since 1900). default [0 inf]
%  dups   0= don't remove any dups
%         1= only remove dups if "primary" datasets are also requested [default]
%         2= remove all "secondary" casts which have dups in "primary" datasets
%            (whether or not the primary datasets are being extracted now.)
%         dups{n} = [m o ...] =remove casts in n duped in datasets m,o,...
%  strp   0= return all casts, even if no data
%         1= return only casts with some data in some profiles [default]
%         2= return only casts with some data in all requested properties
% OUTPUT:
%  lat,lon  [ncast 1]
%  zdep     [ncast ndep]
%  h1 etc   [ncast 1] header vars in order requested
%  v1 etc   [ncast ndep] data vars in order requested
%            {where ndep may be less than number of depths requested, 
%             if no good data that deep for that particular variable.}
%
% Q: Duplicates where more than one variable requested - eg request T & S,
%    one station has T only and the same station in another dataset has S only?   
% A: We have a different list of duplicates for each property, so in the
%    above case both stations would be returned (unless strp=2).
%
% USAGE: [lat,lon,zdep,v1,v2,..] = 
%       get_all_obs(src,hvar,vars,dset,dlim,dmid,nval,scr,trng,dups,strp);

% $Id: get_all_obs.m,v 1.7 2005/11/10 05:44:03 dun216 Exp dun216 $
% Author: Jeff Dunn  CSIRO Marine Research Mar 2005
% Devolved from get_all_csl3.m & get_all_obs.m

function [lat,lon,zdep,varargout] = get_all_obs( ...
    src,hvar,vars,dsets,dlim,dmid,nval,scr,trng,dups,strp)


% MODS:
%  19/6/08  Add access to WOD05
%  22/3/11  Add access to WOD09
%   2/8/11  Add var=9 (pH) extraction for WOD09 (dset 52)
%  25/8/11  Revised TAO format
%  21/3/12  Correcting modified paths

if nargin<4
   error('  GET_ALL_OBS  requires 4 or more input arguments');
end

global MAX_PROF_LEN MAX_PROF_LEN_cnt 

MAX_PROF_LEN = 7000;
MAX_PROF_LEN_cnt = 0;

if nargin<5 
   dlim = [];
elseif length(dlim)==1
   dlim = [0 dlim];
end

if nargin<6
   dmid = [];
end
if nargin<7 | isempty(nval)
   if isempty(dmid)   
      nval = 25000;
   else
      nval = 1;      
   end
end

if nargin<8 | isempty(scr)
   scr = 2;
end
if scr==3
   if (~isempty(dlim) & dlim(1)>0)
      disp('GET_ALL_OBS- WARN: dlim(1)~=0 generally incompatible with scr=3');
   elseif  ~isempty(dmid) 
      disp('GET_ALL_OBS- WARN: use of dmid generally incompatible with scr=3');
   end
end

if nargin<9
   trng = [];
end   

if nargin<10 | isempty(dups)
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

if nargin<11 | isempty(strp)
   strp = 1;
end   

if strp==0 & dups>0
   global Get_All_Obs_told
   if ~Get_All_Obs_told
      disp('You have set dups>0 & strp=0, so duplicate casts will have profiles')
      disp('emptied of data, but they will still be returned');
      Get_All_Obs_told = 1;
   end
end

ii = find(dsets<=6);
if ~isempty(ii)
   disp('This function does not access WOD98');
   dsets(ii) = [];
   if isempty(dsets)
      return
   end
end
   
nhv = length(hvar);
ndv = length(vars);

if min(size(src))==1 & max(size(src))==4
   if src(1)<0; src(1)=0; end
   if src(2)>360; src(2)=360; end
   if src(3)<-90; src(3)=-90; end
   if src(4)>90; src(4)=90; end
end

lat = []; lon = []; zdep = []; cpn = [];

varargout{nhv+ndv} = [];

dupth = path_pc_or_nix('eez_data/boa_duplicates/');

for dset = dsets(:)'
   if any(dset==[21:28 41:49 51:54])
      [lo,la,stnno,zd,iv,hv,vv] = get_wod(dset,src,trng,dlim,nval,dmid,hvar,vars,scr);
   else
      [lo,la,stnno,zd,iv,hv,vv] = get_dset(dset,src,trng,dlim,nval,dmid,hvar,vars,scr);      
   end

   jj = 1:length(lo);
   
   if dups>0 & ~isempty(jj)
      % If dup checking, nan-fill all profiles in duplicate casts.
      % Messy if more than one var, because datasets have different
      % combinations of vars, so vars for the same cast can have
      % different dup status with respect to another dataset.
      % For each var we set dup profiles to nan. We then get rid of whole 
      % stations with only all-nan profiles.
      for jv = iv
	 if dups<3
	    dulst = dset_dup_pref(dupset,vars(jv),dset);
	 else
	    dulst = dupset{dset};
	 end
	 for odset = dulst
	    fnm = sprintf('%s%d_%d_%d_dups',dupth,dset,odset,vars(jv));
	    if exist([fnm '.mat'],'file')
	       load(fnm);
	       ii = inboth(stnno(jj),dupstn);
	       if ~isempty(ii)
		  vv{jv}(jj(ii),:) = nan;
	       end
	    end
	 end	 
      end   
   end


   % If removing casts with all (strp=1) or any (strp=2) empty profiles
   if strp>0 &  ~isempty(jj)
      niv = length(iv);
      if niv==1
	 ii = find(all(isnan(vv{iv}(jj,:)),2));
      else
	 bad = zeros([niv length(jj)]);
	 for jv = 1:niv
	    bad(jv,:) = all(isnan(vv{iv(jv)}(jj,:)),2)';
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
      njj = length(jj);
      ncast = length(lat);

      lat = [lat; la(jj)];
      lon = [lon; lo(jj)];   
      cpn = [cpn; stnno(jj)];
      
      for ii = 1:nhv
	 varargout{ii} = [varargout{ii}; hv{ii}(jj)];
      end

      kk = 1:ndv;
      kk(iv) = [];
      
      % For each variable not available in this dset. Note force minimum 2nd
      % dimension of 1, or will crash when adding a [njj x 0] portion.
      for ii = kk
         ndep = max([1 size(varargout{nhv+ii},2)]);
	 varargout{nhv+ii} = [varargout{nhv+ii}; nan([njj ndep])];
      end
            
      % For each variable in this dset
      for ii = iv
         ndep = size(varargout{nhv+ii},2);
	 ldp = size(vv{ii},2);
	 
	 if ldp < ndep
	    vout = [vv{ii}(jj,:) nan([njj ndep-ldp])];
	    varargout{nhv+ii} = [varargout{nhv+ii}; vout];
	 else
	    if ldp > ndep
	       varargout{nhv+ii} = [varargout{nhv+ii} nan([ncast ldp-ndep])];
	    end
	    varargout{nhv+ii} = [varargout{nhv+ii}; vv{ii}(jj,:)];
	 end
      end
      
      ldp = size(zd,2);
      ndep = size(zdep,2);
      if ldp < ndep
	 zd = [zd(jj,:) nan([njj ndep-ldp])];
	 zdep = [zdep; zd];
      else
	 if ldp > ndep
	    zdep = [zdep nan([ncast ldp-ndep])];
	 end
	 zdep = [zdep; zd(jj,:)];
      end

   end
   
end                                 


if MAX_PROF_LEN_cnt > 0
   disp(['WARNING: get_all_obs truncated ' num2str(MAX_PROF_LEN_cnt) ...
	 ' profiles which exceeded ' num2str(MAX_PROF_LEN) ' values']); 
end

return


%------------------------------------------------------------------------
% Get non-WOD datasets
%
% NOTE: at this stage, only use per-cast QC flags, not the per-value ones.

function [lo,la,stn,zdep,iv,hv,vv] = get_dset(dset,src,trng,dlim,nval,dmid,hvar,vars,scr)

persistent warn13

[pth,slsh] = path_pc_or_nix('eez_data/boa_obslvl/');
pth3c = path_pc_or_nix('eez_data/boa_csl3/csiro/');
infodir = path_pc_or_nix('eez_data/boa_qc_data/');

la = [];
lo = [];
zdep = [];
stn = [];
nhv = length(hvar);
hv{nhv} = [];
vv{length(vars)} = [];

[iv,ih] = dset_vars(dset,1,vars,hvar);
if length(iv)<length(vars)
   if dset==7 && any(vars>2)
      disp('GET_ALL_OBS: dset 7 only has T & S: used dset 72 for nutrients')
   else
      vtmp = vars;
      vtmp(iv) = [];
      disp(['GET_ALL_OBS: vars ' int2str(vtmp) ' not available for dset ' ...
	    int2str(dset)]);      
   end
   if isempty(iv)
      return
   end
end
% Set hvar to 99 where it is not present in this dset, so it will be later
% nan-filled.
hvrs = repmat(99,size(hvar));
hvrs(ih) = hvar(ih);

vnm = {'t','s','o2','si','po4','no3','nutdens','','','','','','','no4','nh3'};

deps = [];
stdmf = 0;
stdnc = 0;

switch dset
           
  case 7 
    % CSIRO obs level data is stored in:
    % - CTD-only CSLv3 file - accessed here, but only has T&S
    % - 2db per-cruise files - see dset 70
    % - hydro obs file - see dset 72
    % - historical no3 obs file - see dset 71
    % Combined CTD and hydro obs CSLv3 file "csiro_csl.nc" is used for CSL3 
    % access only.
    
    cellstored = 0;
    stdnc = 1;    
    fnm = [pth3c 'csiro_ctd_csl'];
    lat = getnc(fnm,'lat');
    lon = getnc(fnm,'lon');
    time = getnc(fnm,'time');
    stnno = getnc(fnm,'stnno');
    cdep = getnc(fnm,'depth');
    if any(hvrs==4)
       botdep = getnc(fnm,'botdepth');
    end
    if any(hvrs==5)
       co = repmat(12345,[length(lat) 1]);
    end
    if any(hvrs==6)
       cru = floor(stnno/1000);
    end
      
            
  case 9
    % NIWA
    %cellstored = 0;
    %stdnc = 1;
    %fnm = [pth 'niwa' slsh 'niwa_obs'];
    %lat = getnc(fnm,'lat');
    %lon = getnc(fnm,'lon');
    %time = getnc(fnm,'time');
    %stnno = getnc(fnm,'stnno');
    %cru = floor(stnno/10000);
    %pp = getnc(fnm,'pressure');
    %deps = nan(size(pp));
    %co = repmat(13873,size(lat));
    %for jj = 1:length(lat)
    %   deps(jj,:) = sw_dpth(pp(jj,:),lat(jj));
    %end
    
    cellstored = 1;
    stdmf = 1;
    load([pth 'niwa' slsh 'niwa_obs']);
    co = repmat(13873,size(lat));
        
  case 10 
    % French Indian 2001 CD
    cellstored = 1;
    stdmf = 1;
    load([pth 'french_data' slsh 'french_hydro']);
    deps = depth';
    botdep = botdepth;
    % Cruise codes are crap, but this is the best we have..
    cru = cru_num(:);

  
  case 11
    % Argo
    cellstored = 0;
    stdmf = 1;
    load([pth 'argo' slsh 'argo_obs']);
    
  case {12,61,62}
    % TAO
    cellstored = 0;
    stdmf = 1;
    if dset==12
       LL = load([pth 'tao' slsh 'tao_obs']);
    elseif dset==61
       LL = load([pth 'tao' slsh 'pirata_obs']);
    elseif dset==62
       LL = load([pth 'tao' slsh 'rama_obs']);
    end
    if any(vars==1) & any(vars==2)
       disp('** TAO has different depths for T and S - so each variable')
       disp('must be extracted separately. ***');
    end
    if any(vars==1)
       cdep = LL.t_depth;
       sz = size(LL.t);
       stnno = vec(LL.stnno');
       tim = LL.time;
       t = zeros(prod(sz([1 3])),sz(2));
    else
       cdep = LL.s_depth;
       sz = size(LL.s);
       tim = LL.time(LL.itim_t_s);
       stnno = vec(LL.stnno(LL.itim_t_s,:)');
       s = zeros(prod(sz([1 3])),sz(2));
    end
    ntot = prod(sz([1 3]));
    lon = zeros(ntot,1);
    lat = zeros(ntot,1);
    cru = zeros(ntot,1);
    time = zeros(ntot,1);
    botdep = zeros(ntot,1);
    for ij = 1:sz(1)
       jk = ((ij-1)*sz(3)) + (1:sz(3));
       lon(jk) = LL.lon(:);
       lat(jk) = LL.lat(:);
       if isfield(LL,'cru')
	  cru(jk) = LL.cru(:);
       else
	  cru(jk) = 1:sz(3);
       end
       botdep(jk) = LL.botdep(:);
       time(jk) = tim(ij);
       if any(vars==1)
	  t(jk,:) = sq(LL.t(ij,:,:))';
       else
	  s(jk,:) = sq(LL.s(ij,:,:))';
       end
    end
    clear LL tim
  
  case 13
    % Antarctic CRC - 
    % This file is the worst. T,S,o2 and nutdens are in CSLv3 matrices; no3,
    % si and po4 are in cells and their z is pressure in var obspres. If want
    % some of both, have to convert all to cells, on the fly, combining the 
    % CSL and obs depths. This also means that per-value screening indices 
    % into profiles will be mucked up if access both at once. 
    
    load([pth 'crc' slsh 'crc_csl_obs']);    
    co = repmat(12345,size(stnno));
    cru = cru_idx;
    
    ivm = find(ismember(vars,[1 2 3 7]));
    ivc = find(ismember(vars,[4 5 6]));
    cdep = csl_dep(1:79,3);    

    cellstored = ~isempty(ivc);
    
    if cellstored
       if scr>1 && ~isempty(ivm) && isempty(warn13) 
	  warn13 = 1;
	  disp(['** WARNING: scr>=2 conflicts with this combination of ' ...
		'properties for dset 13']);
       end
       clear deps
       for ii = 1:length(lat)
	  if isempty(obspres{ii})
	     deps{ii} = cdep;
	     for jj = ivm
		eval(['vin{ii}{jj} = ' vnm{vars(jj)} '(ii,:);']);
	     end
	  else
	     obsdep = sw_dpth(obspres{ii},lat(ii));
	     [deps{ii},indx] = sort([cdep(:); obsdep(:)]);
	     blnk = nan([1 length(deps{ii})]);
	     kk = find(indx>79);
	     for jj = ivc
		eval(['vtmp = ' vnm{vars(jj)} '{ii};']);
		if ~isempty(vtmp)
		   vin{jj}{ii} = blnk;
		   vin{jj}{ii}(kk) = vtmp;
		end
	     end
	     kk = find(indx<=79);
	     for jj = ivm	     
		vin{jj}{ii} = blnk;
		eval(['vin{jj}{ii}(kk) = ' vnm{vars(jj)} '(ii,:);']);
	     end
	  end
       end
       % Now ensure that we have at least an empty cell for every station 
       % and property
       for jj = 1:length(vin)
	  if length(vin{jj})<length(lat)
	     vin{jj}{length(lat)} = [];
	  end
       end
    else       
       for ii = ivm
	  eval(['vin{ii} = ' vnm{vars(ii)} ';']);  
       end
    end
    for ii = iv	  
       if vars(ii)==7
	  cflg{ii} = zeros(1,length(lat));
       else
	  eval(['cflg{ii} = ' vnm{vars(ii)} '_castflag;']);
       end
    end       
  
  
  case 14
    % Willis global
    % This is a huge and homogeneous dataset (all casts to 750m, no nans)
    % Easy to call get_willis directly, and that way uses less memory! 
    cdep = 0:10:750;
    idp = req_depth(dlim,nval,dmid,cdep);
    if scr==3 & idp(1)~=1
       disp('GET_ALL_OBS - WARNING scr=3 incompatible with dlim for dset=14'); 
    end
    cdep = cdep(idp);
    [lon,lat,time,t,botdep,a1,a2,a3,a4,stnno] = get_willis(src,cdep);
    clear a?
    
    
  case 19
    % WOCE3 WHP - this is pretty complicated. Have file WHPheads.mat, then
    % woce3_ctd_ p, deps, o, s, t, then woce3_hyd. woce3_hyd contains bot_s,t 
    % as well as hyd_s,t, also vars nite & nate.  For now, coded access to CTD
    % data; leave hydro for now in the hope I will not have to code it!
    cellstored = 1;
    stdmf = 1;
    load([pth 'woce_whp' slsh 'woce3_ctd_deps']);
    load([pth 'woce_whp' slsh 'WHPheads'],'stnno','time','cri1','lat','lon');
    cru = cri1;
    if any(hvrs==4)
       load([pth 'woce_whp' slsh 'woce3_ctd_p'],'botdep');
    end 
    if any(vars==1)
       load([pth 'woce_whp' slsh 'woce3_ctd_t']);
    end
    if any(vars==2)
       load([pth 'woce_whp' slsh 'woce3_ctd_s']);
    end
    if any(vars==3)
       load([pth 'woce_whp' slsh 'woce3_ctd_o2']);
    end
    if any(vars>3)
       disp('GET_ALL_OBS: no access yet coded to WOCE hydro variables')
    end
    
    
  case 31
    % Other - mainly AIMS Torres CTD 
    load([pth 'other' slsh 'other_obslvl']);
    cellstored = 1;
    stdmf = 1;    
    
  case 35
    % AWI (Alfred Wegener Institute) CTD
    cellstored = 1;
    stdmf = 1;
    load([pth 'awi' slsh 'awi_obs']);
    botdep = bdep;
    
  case 70
    % CSIRO 2db CTD - per-cruise files. Could be accessed here, but will 
    % not code that unless shown to be required.
    disp('CSIRO 2db CTD per-cruise files are not presently accessible via this')
    disp('function.  Use read_csiro_2db.m instead.')
    
    
  case 71
    % CSIRO historical nitrate
    cellstored = 1;
    stdmf = 1;
    load([pth 'csiro' slsh 'csiro_no3_obs']);
    stnno = prid;
    if any(hvrs==4)
       botdep = bdpth;
    end
    if any(hvrs==6)
       cru = floor(stnno/1000);
    end

    
  case 72
    % CSIRO Hydro Obs-level file    
    cellstored = 1;
    stdmf = 1;
    load([pth 'csiro' slsh 'csiro_hyd_obs']);
    botdep = bdpth;
    if any(hvrs==5)
       co = repmat(12345,[length(lat) 1]);
    end
    if any(hvrs==6)
       cru = floor(stnno/1000);
    end
    
end


if stdnc
   for ii = iv
      vin{ii} = getnc(fnm,vnm{vars(ii)});
      if vars(ii)==7
	 cflg{ii} = zeros(1,length(lat));
      else	 
	 cflg{ii} = getnc(fnm,[vnm{vars(ii)} '_castflag']);
      end
   end  
end

if stdmf
   for ii = iv
      eval(['vin{ii} = ' vnm{vars(ii)} '; clear ' vnm{vars(ii)} ';']);      
      if exist([vnm{vars(ii)} '_castflag'],'var')
	 eval(['cflg{ii} = ' vnm{vars(ii)} '_castflag;']);
      else
	 cflg{ii} = zeros(1,length(lat));
      end
   end
   clear *_castflag
end


if isempty(src)
   jj = 1:length(lon);
elseif min(size(src))==1
   jj = find(lon>=src(1) & lon<=src(2) & lat>=src(3) & lat<=src(4));
else
   jj = find(inpolygon(lon,lat,src(:,1),src(:,2)));
end
    
if ~isempty(trng) & ~isempty(jj)
   kk = find(time(jj)<trng(1) | time(jj)>=trng(2));
   jj(kk) = [];
end


if ~isempty(jj)
   
   % Initialising variables
   ntot = 0;
   nblk = 500;
   zdep = nan([nblk 1]);
   for ii = iv
      vv{ii} = nan([nblk 1]);
   end
   nbuf = nblk;
    
   % Before finding and loading data, apply the required level of screening
   % by nan-ing any flagged data (or removing that data from cells). 
   % Pre-working this saves dealing with lots of cases while finding and
   % loading the data.  
   
   for ii = iv
      % For each property, per-value screen if required:
      if scr>=2 & vars(ii)<=6 & ~isempty(vv{ii}) 
	 pnm = prop_name(vars(ii));
	 fnm = [infodir num2str(dset) slsh pnm '_obs_scr'];
	 if exist([fnm '.mat'],'file')
	    load(fnm);
	    [tmp,ijj,iscr] = intersect(stnno(jj),scrstn);
	    
	    % Tribute to the god of indexing....
	    if cellstored & ~isempty(ijj)
	       for kk = 1:length(ijj)
		  vin{ii}{jj(ijj(kk))}(scridx{iscr(kk)}) = nan;
	       end
	    elseif  ~isempty(ijj)
	       for kk = 1:length(ijj)
		  vin{ii}(jj(ijj(kk)),scridx{iscr(kk)}) = nan;
	       end
	    end		  
	 end
      end
	 
      % Per-cast screening, if requested
      if scr>=1 
	 kk = jj(find(cflg{ii}(jj)~=0));
	 if cellstored
	    for jk = kk(:)'
	       vin{ii}{jk} = [];
	    end	    
	 else
	    vin{ii}(kk,:) = nan;
	 end
      end            
   end  % Loop on iv


   if cellstored   
      for nn = 1:length(jj)
	 kk = jj(nn);
	 % For each required station, find all depths in required range which
	 % have valid data for any of the required properties.
	 
	 ll = [];
	 for ii = iv
	    if ~isempty(vin{ii}{kk})
	       if isempty(ll)		  
		  ll = ~isnan(vin{ii}{kk});
	       else
		  ll = ll | ~isnan(vin{ii}{kk});
	       end
	    end
	 end
	 
	 if any(ll)
	    if scr==3
	       ll = 1:length(deps{kk});
	    else
	       ll = find(ll);
	    end
	    idp = ll(req_depth(dlim,nval,dmid,deps{kk}(ll)));
	 else
	    idp = [];
	 end	 

	 % We will reject all casts which have no data in the required depth range 
         if isempty(idp)
	    jj(nn) = 0;
	 else
	    ntot = ntot+1;
	    if ntot>nbuf
	       % then need to extend pre-allocation of output arrays
	       for ii = iv
		  vv{ii} = [vv{ii}; nan([nblk size(vv{ii},2)])];
	       end
	       zdep = [zdep; nan([nblk size(zdep,2)])];
	       nbuf = nbuf+nblk;
	    end
	    
	    ndep = length(idp);	    
	    if ndep > size(zdep,2)
	       zdep = [zdep nan([nbuf ndep-size(zdep,2)])];
	    end
	    zdep(ntot,1:ndep) = deps{kk}(idp);
	       
	    for ii = iv
	       if ~isempty(vin{ii}{kk})
		  if ndep>size(vv{ii},2)
		     vv{ii} = [vv{ii} nan([nbuf ndep-size(vv{ii},2)])];
		  end
		  vv{ii}(ntot,1:ndep) = vin{ii}{kk}(idp);
	       end
	    end
	 end
      end

   else
      % Stored as matrices
      
      % For now, check vars are dimensioned properly
      for ii = iv
	 if size(vin{ii},1) ~= length(lon)
	    disp(['Transpose var ' num2str(vars(ii)) ' in dset ' num2str(dset)]);
	    vin{ii} = vin{ii}.';
	 end
      end

      if isempty(deps)
	 deps = repmat(cdep(:)',[max(jj) 1]);
      end
      ldep = size(deps,2);
      
      for nn = 1:length(jj)
	 kk = jj(nn);
	 % For each required station, find all depths in required range which
	 % have valid data for any of the required properties.

	 ll = zeros(1,ldep);
	 for ii = iv
	    if any(~isnan(vin{ii}(kk,:)))
	       ll = ll | ~isnan(vin{ii}(kk,:));
	    end
	 end

	 if any(ll)
	    if scr==3
	       ll = 1:ldep;
	    else
	       ll = find(ll);
	    end
	    idp = ll(req_depth(dlim,nval,dmid,deps(kk,ll)));
	 else
	    idp = [];
	 end	 

         if isempty(idp)
	    % will reject casts which have no data in the required depth range
	    jj(nn) = 0;
	 else
	    ntot = ntot+1;

	    if ntot>nbuf
	       % then need to extend pre-allocation of output arrays
	       for ii = iv
		  vv{ii} = [vv{ii}; nan([nblk size(vv{ii},2)])];
	       end
	       zdep = [zdep; nan([nblk size(zdep,2)])];
	       nbuf = nbuf+nblk;
	    end
	    
	    ndep = length(idp);
	    if ndep > size(zdep,2)
	       zdep = [zdep nan([nbuf ndep-size(zdep,2)])];
	    end
	    zdep(ntot,1:ndep) = deps(kk,idp);
	    for ii = iv
	       if ndep>size(vv{ii},2)
		  vv{ii} = [vv{ii} nan([nbuf ndep-size(vv{ii},2)])];
	       end
	       vv{ii}(ntot,1:ndep) = vin{ii}(kk,idp);
	    end	 
	 end
      end
   end   
   
   jj(find(jj==0)) = [];   
   nj = length(jj);

   if nj>0
      la = lat(jj);
      lo = lon(jj);
      stn = stnno(jj);
   
      for ii = 1:nhv
	 switch hvrs(ii)
	   case 1
	     hv{ii} = stn;
	   case 2
	     hv{ii} = time(jj);
	   case 3
	     hv{ii} = cflg{1}(jj);
	   case 4
	     hv{ii} = botdep(jj);
	   case 5
	     hv{ii} = co(jj);
	   case 6
	     hv{ii} = cru(jj);
	   case 7
	     if exist('dmode','var')
		hv{ii} = dmode(jj);
	     end
	   case 9
	     hv{ii} = repmat(dset,[nj 1]);
	   case 13
	     if exist('rejpc','var')
		hv{ii} = rejpc(jj);
	     end
	   otherwise
	     hv{ii} = nan([nj 1]);		 
	 end
      end	
   end   

   % Trim pre-allocated variables
   zdep = zdep(1:ntot,:);
   for ii = iv
      vv{ii} = vv{ii}(1:ntot,:);
   end
end


return

%------------------------------------------------------------------------
% WOD is different because stored in WMO files.

function [lo,la,pn,zdep,iv,hv,vv] = get_wod(dset,src,trng,dlim,nval,dmid,hvar,vars,scr);

global MAX_PROF_LEN MAX_PROF_LEN_cnt 

la = [];
lo = [];
pn = [];
zdep = [];
nhv = length(hvar);
hv{nhv} = [];
vv{length(vars)} = [];

if any(dset==[25 26 27 28 29 45 47 48 49 55:59])
   disp(['No access to dataset ' num2str(dset) ' yet!'])
   return
end
  
iv = dset_vars(dset,1,vars);

% Conversion between our property codes and these used in the WOD files:
%        t  s  o2  si  po4  no3  -  -  pH 
wvar =  [1  2  3   6    4    8   0  0  9];

wod05 = (dset>40 & dset<50);
wod09 = dset>50;

[pth,slsh] = path_pc_or_nix('eez_data/boa_obslvl/wod0');
if wod09
   [pth,slsh] = path_pc_or_nix('datalib/observations/WOD09/obs_level/');
   prefx = {'CTD','OSD','PFL','UOR','DRB','GLD','MRB','XBT'};
   ids = dset-50;
elseif wod05
   pth = [pth '5' slsh];
   prefx = {'CTD','OSD','PFL','UOR','DRB','GLD','MRB','XBT'};
   ids = dset-40;
else
   pth = [pth '1' slsh];
   prefx = {'CTD','OSD','PFL','UOR','DRB','MBT','MRB','XBT'};
   ids = dset-20;
end

if isempty(src)
   src = [0 360 -80 80];
end    
if min(size(src))==1
   wmosq = getwmo(src);
else
   wmosq = getwmo(...
       [min(src(:,1)) max(src(:,1)) min(src(:,2)) max(src(:,2))]);
end


% Initialise variables
ntot = 0;
nblk = 500;
zdep = nan([nblk 1]);
for ii = iv
   vv{ii} = nan([nblk 1]);
end
nbuf = nblk;


for wmo = wmosq
   flnm = [pth prefx{ids} slsh prefx{ids} num2str(wmo)];

   if ~exist([flnm '.mat'],'file')
      jj = [];
   else
      load(flnm)
      if min(size(src))==1
	 jj = find(lon>=src(1) & lon<=src(2) & lat>=src(3) & lat<=src(4));
      else
	 jj = find(inpolygon(lon,lat,src(:,1),src(:,2)))';
      end

      if ~isempty(trng) & ~isempty(jj)
	 kk = find(time(jj)<trng(1) | time(jj)>=trng(2));
	 jj(kk) = [];
      end
   end
   
   if ~isempty(jj)      
      % If screening, clobber property codes whereever "cast flag" is set
      ipsave = ipcode;
      if scr>=0
	 rr = find(ierror~=0);
	 ipcode(rr) = 0;
      end
      
      for nn = 1:length(jj)
	 kk = jj(nn);
	 % For each required station, find all depths in required range which
	 % have valid data for any of the required properties. 
	 % ivw indexes vars which are present for this cast
	 % wiv indexes the profs which match those vars 
	 ll = zeros(1,length(deps{kk}));
	 ivw = [];
	 for ii = iv
	    gg = find(ipcode(kk,:)==wvar(vars(ii)));
	    %## if ~isempty(gg)
	    if ~isempty(gg) && size(profs{kk},2)<gg
	       disp(num2str([wmo size(profs{kk},2) gg]))
	    elseif ~isempty(gg)
	       if scr>=2 
		  ll = ll | (~isnan(profs{kk}(:,gg)) & iderror{kk}(:,gg)==0)';
	       else	
		  ll = ll | ~isnan(profs{kk}(:,gg))';
	       end
	       ivw = [ivw ii];
	       wiv(ii) = gg;
	    end
	 end

	 if any(ll)
	    if scr==3
	       ll = 1:length(deps{kk});
	    else
	       ll = find(ll);
	    end
            idp = ll(req_depth(dlim,nval,dmid,deps{kk}(ll)));
	    if length(idp) > MAX_PROF_LEN
	       idp = idp(1:MAX_PROF_LEN);
	       MAX_PROF_LEN_cnt = MAX_PROF_LEN_cnt + 1;
	    end
	 else
            idp = [];
         end

         if isempty(idp)
	    % will reject casts which have no data in the required depth range
	    jj(nn) = 0;
	 else
            ntot = ntot+1;
            if ntot>nbuf
               for ii = iv
                  vv{ii} = [vv{ii}; nan([nblk size(vv{ii},2)])];
               end
               zdep = [zdep; nan([nblk size(zdep,2)])];
               nbuf = nbuf+nblk;
            end

	    ndep = length(idp);
	    if ndep > size(zdep,2)
	       zdep = [zdep nan([nbuf ndep-size(zdep,2)])];
	    end
            zdep(ntot,1:ndep) = deps{kk}(idp);

            for ii = ivw
               if ndep>size(vv{ii},2)
                  vv{ii} = [vv{ii} nan([nbuf ndep-size(vv{ii},2)])];
               end
	       if scr>=2 & any(iderror{kk}(idp,wiv(ii)) > 0)
		  mm = find(iderror{kk}(:,wiv(ii)) > 0);
		  profs{kk}(mm,wiv(ii)) = nan;
	       end
	       vv{ii}(ntot,1:ndep) = profs{kk}(idp,wiv(ii))';
            end
         end
      end       % for each suitable cast

      jj(find(jj==0)) = [];
   end	% if some suitable casts in this file 
	 
	 
   if ~isempty(jj)
      nj = length(jj);
	 
      if size(lat,1)==1
	 lat = lat';
	 lon = lon';
	 ocl = ocl';
	 lat = lat';
	 time = time';
	 botdep = botdep';
	 cc = cc';
	 cru = cru';
      end
	 
      la = [la; lat(jj)];
      lo = [lo; lon(jj)];
      pn = [pn; ocl(jj)];
	 
      for ii = 1:nhv
	 switch hvar(ii)
	   case 1
	     hv{ii} = [hv{ii}; ocl(jj)];
	   case 2
	     hv{ii} = [hv{ii}; time(jj)];
	   case 3
	     tmp = zeros([nj 1]);
	     for kj = 1:length(jj)
		lj = find(ipsave(jj(kj),:)==wvar(vars(1)));
		if ~isempty(lj)
		   tmp(kj) = ierror(jj(kj),lj);
		end
	     end
	     hv{ii} = [hv{ii}; tmp];
	   case 4
	     hv{ii} = [hv{ii}; botdep(jj)];
	   case 5
	     hv{ii} = [hv{ii}; cc(jj)];
	   case 6
	     hv{ii} = [hv{ii}; cru(jj)];
	   case 9
	     hv{ii} = [hv{ii}; repmat(dset,[nj 1])];
	   otherwise
	     hv{ii} = [hv{ii}; nan([nj 1])];		 
	 end
      end	      
   end    % endif ~isempty(jj)

end     % endfor wmo = wmosq 


% Trim pre-allocated variables
zdep = zdep(1:ntot,:);
for ii = iv
   vv{ii} = vv{ii}(1:ntot,:);
end



%-----------------------------------------------------------------------------
function idp = req_depth(dlim,nval,dmid,deps)

            
% Limit to depth range, if required
if ~isempty(dlim)
   idp = find(deps>=dlim(1) & deps<=dlim(2));
else
   idp = 1:length(deps);
end

% Limit number of values, if required
if length(idp)>nval
   % Extract around a mid-depth, if required
   if ~isempty(dmid)
      [tmp,sss] = sort(abs(deps(idp)-dmid));
      idp = idp(sort(sss(1:nval)));
   else
      idp = idp(1:nval);
   end
end

%---------------------------------------------------------------------------
