% GET_SST_XY  Get SST at given locations and for a *SINGLE TIME*, from any of 
%    the local datasets. The nearest time estimate is returned, except for 
%    datasets 22, 28, 29, 30, 31, 33 where we interpolate to the specified time.
%
% SEE http://www.marine.csiro.au/eez_data/doc/sst_datasets.html
%
% In early 2007 that whole /home/satdata1/ directory was abolished, and the 
% "Decade of SST" ("stitiched") files were setup for permanant updating.
%
% Note that this code is very slow - it spends most time working out which 
% file to open next, so especially if a time series is required, manually
% accessing the 3D netcdf files is much quicker. 
%
% INPUT
%  x,y  - vector or matrices of locations
%  tim  - single Unix time of required data; eg greg2time([2001 12 31 0 0 0])
%  pref - vector of dataset preferences (default [18 22 35 16 17])
%        Note that unless option 1 is used, all datasets will be tried until
%        there are no missing values, but data is always missing over land, 
%        so all pref datasets *will* be slowly accessed if there is land
%        at any x,y.
%        Data prior to 1/10/93 is available only in Reynolds analyses, but 
%        these are not default datasets because they are a different type
%        of product (but by all means use them if it suits you.)
%         1  Walker/Wilkin Pathfinder space-time 9km 10day OI
%         2  [ DEFUNCT ]  Walker/Wilkin ACRES Pathfinder space-time OI 5km 10day OI
%         3  [ DEFUNCT ]  3-day comp east
%         5  [ DEFUNCT ]  5-day comp west
%         6  [ DEFUNCT ]  6-day comp east
%         7  [ DEFUNCT ]  10-day comp southern 
%         8  [ DEFUNCT ]  Rathbone/Griffin east composite
%         9  [ DEFUNCT ]  Rathbone/Griffin west composite
%        10  [ DEFUNCT ]  10-day comp east
%        11  15-day comp GAB - RESTORED 30/12/1989 to 31/12/1994
%        12  [ DEFUNCT ]  10-day comp west 
%        13  [ DEFUNCT ]  10-day comp full Aus  20/7/02 ->
%        14  [ DEFUNCT ]  6-day comp  full Aus  11/7/02 ->
%        15  [ DEFUNCT ]  3-day comp  full Aus  19/7/02 ->
%        16  3-day ave AMSR SSMI Microwave, Global 1/6/02 -> 4/10/2011
%        17  7-day ave AMSR SSMI Microwave, Global 1/6/02 -> 1/10/2011
%   SST AVHRR Australasian 4 km resolution time series
%   (previously called "2004 AVHRR NOO Stitched")
%                                  Yearly and 6-monthly files
%        18  1-day  1-day spaced 1/10/93->
%        19  3-day  1-day spaced 1/10/93->
%        20  6-day  1-day spaced 1/10/93->
%        21  10-day 2-day spaced 1/10/93->
%        22  15-day 6-day spaced 1/10/93->
%            Daily files (previously called "AVHRR Aasia NRT")
%            [these are actually available to earlier dates but in a
%            different location - talk to Jeff if they are required.]
%        23  3-day  20/10/04 ->  
%        24  6-day  19/10/04 ->
%        25  10-day 17/10/04 ->
%        26  15-day   [not available?]
%        27  20-day   [not available?]
%        28  6-day "Patchwork" V7  1/6th degree 15/6/92-17/8/2004
%        29  6-day "Patchwork" V8  1/6th degree 17/8/2004-26/4/2005
%        30  6-day Reynolds V7  1 degree 15/6/92-17/8/2004
%        31  6-day Reynolds V8  1 degree 17/8/2004-26/4/2005
%        33  7-day Reynolds V2  1 degree 29/10/1981 ->
%        34  Monthly NOAA V3 2 degree 1854 ->
% *NEW*  35  Windsat v7 daily  0.25 degree  1/1/2011 ->   (may take over from 16,17)
%            [See opt=5, used to specify descending or ascending passes, or both]
%
%  opt  - vector of any of the following options (defaults in []):
%         1: 1=get from one dataset only, not top-up from others [0]
%         2: 1=single timu value from first dataset used [0]
%         3: min "nclear" required for post-2004 AVHRR 4km yearly file data
%             A rough guide from David Griffin (for 20 day composites?):
%            5  - for risk-takers
%            10 - for prudent oceanographers
%            20 - for perfectionists
%            [Default: no test of nclear]  
% *NEW*   4: 1=spread edges of data before interpolation (because if even one
%            surrounding point is NaN, interp result will be NaN, so high losses if
%            scattered bad points. This option takes much longer! [0]  
%         5: 0=all passes  1=ascending only 2=descending only (dset 35 only) [0]
%         
%  o_val - vector of values corresponding to the options specified above
%
% OUTPUT
%  sst  - SST at x,y locations (nan where no data available)
%  timu - timestamp of data used (nth element is value for nth dataset, as
%         listed above - eg if only data from "3-day comp", the only non-zero
%         value will be in element 3). Note that this refers to centre of
%         time window.
%  dset - dataset used (see pref above) for each x,y
%  
% SEE ALSO:   get_sst_xyt.m   get_sst.m
%
% CALLS:      sst_filnam.m
%
% Jeff Dunn  CSIRO CMR 5/4/02 ... 23/3/10
% 
% USAGE: [sst,timu,dset] = get_sst_xy(x,y,tim,pref,opt,o_val);

function [sst,timu,dset] = get_sst_xy(x,y,tim,pref,opt,o_val)

% MODS:    Now access netCDF rather than disimp daily files       28/2/03
%          Access to ACRES Pathfinder 5km 10day OI                4/3/03
%          Added Global 3-day & 7-day SSMI (AMSR) microwave SST    8/8/03
%          Access to 2004 Stitched Archive composites             10/6/04
%          Access to AVHRR NRT composites                         18/7/05
%          Patchwork V8 and Reynolds V7,V8.                       1/8/06
%          Reynolds V1,V2, remove refs to deleted datasets, 
%          handle change from yearly to 6monthly NOO files        14/8/07
%          Testing of nclear variable for AVHRR 4km yearly files  28/8/07
%          Restored lost dataset 11                               8/4/08
%          Moving 30-33, add 34                                   23/3/10
%          Add 35, and edge-spread option                         23/2/12
%          Add option 5                                           23/3/12
%          Check/update file locations                            7/11/13
%
% ON DATASETS...
% A reorganisation around 20/7/02 lead to many datasets being terminated,
% and the full Aus Region ("Au") ones commencing at 3, 6, and 10 day windows
% in /home/satdata1/SST/.  All composite data was then available in disimp,
% single-date netcdf, and 3D netcdf files. In May 2004 the NOO-funded
% "Decade of SST" project produced 1,2,3,10 and 15-day window AVHRR composites 
% of the Australasian region in netCDF year files. However, many of the
% previous datasets have been retained as they extend beyond the Decade, or 
% have unique characteristics. 
%
% In early 2007 that whole /home/satdata1/ directory was abolished, and the 
% "Decade of SST" ("stitiched") files were setup for permanant updating.

sst=[]; timu=[]; dset = [];

% Number of underscores in filename before start of date string (for each
% data stream)
ndsh = [1 1 1 1 1 1 1 1 1 1 0 1 1 1 1 3 2 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0];

% Days to subtract from file date to get analysis or window centre date
% (for each data stream.)  [So far, have only considered this for 16 & 17]
toff = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 .5 2.5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
	0 0];

if nargin<4 || isempty(pref)
   pref = [18 22 35 16 17];
elseif any(ismember(pref,[2:10 12:15]))
   kk = find(ismember(pref,[2:10 12:15]));
   disp(['** GET_SST_XY no longer accesses datasets ' num2str(pref(kk))]);
   pref(kk) = [];
elseif any(pref==32)
   kk = find(pref==32);
   disp(['GET_SST_XY: use dset 33 instead of 32, or directly access the' ...
	 ' monthly Reynolds V2 file.']);
   pref(kk) = [];
elseif any(pref==34)
   kk = find(pref==34);
   disp(['GET_SST_XY: Contact Jeff Dunn to code access to Dset 34. I do not' ...
	 ' intend to code it unless someone actually needs it!']);
end
if isempty(pref)
   return
end
npref = length(pref);

% What sort of computer are we running on?
cname = computer;
if strncmp(cname,'PC',2)
   pth = '\\rsl\rsl3\datasets\avhrr\sstcr_old\SOUTH_oz\netcdf\';
   pth2 = '\\reg2\reg2\SST_mw\netcdf\';
   pth3 = '\\rsm\imgjj\sensor\avhrr\sstcr04\yearfiles\';
   pth4 = '\\rsm\imgjj\sensor\avhrr\sstcr04nrt\';
   pth4b = '\\rsm\imgjj\sensor\avhrr\sstcr04\';
   pth5 = '\\fstas2-hba\datalib\reanalyses\Reynolds\';   
   pth6 = '\\fstas2-hba\CMAR-HOME1\pigot\sst_nc_wm\';   
   pth7 = '\\fstas2-hba\datalib\reanalyses\Reynolds\Reconstructed\';   
   pth8 = '\\fstas2-hba\datalib\platforms\wsat\unzipped\';
   slsh = '\';
elseif strncmp(cname,'MAC',3)
   disp([7 'Sorry - do not how to find datafiles from a Mac'])
   return
else
   pth = '/home/rsl3/datasets/avhrr/sstcr_old/SOUTH_oz/netcdf/';
   pth2 = '/home/reg2/SST_mw/netcdf/';
   pth3 = '/home/imgjj/sensor/avhrr/sstcr04/yearfiles/';
   pth4 = '/home/imgjj/sensor/avhrr/sstcr04nrt/';
   pth4b = '/home/imgjj/sensor/avhrr/sstcr04/';
   pth5 = '/home/datalib/reanalyses/Reynolds/';
   pth6 = '/home/pigot/sst_nc_wm/';
   pth7 = '/home/datalib/reanalyses/Reynolds/Reconstructed/';
   pth8 = '/home/datalib/platforms/wsat/unzipped/';
   slsh = '/';
end
   
sngl = 0;
tmu1 = 0;
mincl = [];
sprd = 0;
passprf = 0;
if nargin<5 | isempty(opt)
   % use defaults above
elseif nargin<6 | isempty(o_val)
   warning('GET_SST_XY: Need to specify values to go with "opt" options');
else
   for ii = 1:length(opt)
      switch opt(ii)
	case 1
	  sngl = o_val(ii);
	case 2
	  tmu1 = o_val(ii);
	case 3
	  mincl = o_val(ii);
	case 4
	  sprd = o_val(ii);
	case 5
	  passprf = o_val(ii);
	otherwise
	  disp(['GET_SST_XY: Option ' num2str(opt(ii)) ' means nix to me!']);
      end
   end
end

% Create output variables and set up loop control
dset = zeros(size(x));
timu = zeros([1 max(pref)]);
sst = repmat(nan,size(x));
nval = prod(size(x));
jj = 1:nval;
req = 1;
ip = 1;
nclr = [];

% Loop through the possible datasets in order of preference

while req
   pr = pref(ip);
   
   t0 = greg2time([1970 1 1 0 0 0]);
   
   % Get time and region limits for this dataset 
   switch pr
     case 1
       % Walker/Wilkin PF 9km 10-day OI  14/2/1989-27/6/1994
       t1 = 31820; 
       t2 = 34510;
       gg = [90 199.95 -69.96 -.09];
     %case 2
       % Walker/Wilkin PF 10-day OI ACRES 5km 18/3/1991-3/8/1997
      % t1 = 33313;
      % t2 = 35643;
      % gg = [106.80 160.16 -44.94 -3.95];
      % fnm = [pth 'oi5k10dsst1'];
      % t0 = greg2time([1990 1 1 0 0 0]);
     %case 8
       % Rathbone/Griffin east - 15/9/1991-13/3/1999
      % t1 = 33484;
      % t2 = 36230;
      % gg = [147.26 156.92 -43.02 -25.32];
      % fnm = [pth 'comp15deast'];
      % t0 = greg2time([1990 1 1 0 0 0]);
     % OLD case 11
       % OLD case 11: 15-day comp GAB - 7/12/1989-17/7/02
      % t1 = 32847;
      % t2 = 37452;
      % gg = [110.01 159.96 -46.96 -30.01];
      % dnm = [pth 'comp15d' slsh 'GAB' slsh 'nc' slsh];
     case 11
       % Restored 15-day comp GAB, 5 day spacing - 30/12/1989-23/3/2004
       % Filename date refers to middle of time window.
       t1 = 32870;
       t2 = 38067;
       gg = [110.01 159.97 -46.97 -30.01];
       dnm = pth;
     %case 13
       % 10-day comp Aus - 20/7/02 - present
      % t1 = 37455;
      % t2 = inf;
      % gg = [100.00 171.42 -47.03 -0.018];
      % dnm = [pth 'comp10d' slsh 'Au' slsh 'nc' slsh];
     %case 14
       % 6-day comp Aus - 11/7/02 - present
      % t1 = 37446;
      % t2 = inf;
      % gg = [100.00 171.42 -47.03 -0.018];
      % dnm = [pth 'comp6d' slsh 'Au' slsh 'nc' slsh];
     %case 15
       % 3-day comp Aus - 19/7/02 - present
      % t1 = 37454;
      % t2 = inf;
      % gg = [100.00 171.42 -47.03 -0.018];
      % dnm = [pth 'comp3d' slsh 'Au' slsh 'nc' slsh];
      
     case 16
       % 3-day ave SSMI Microwave;  Global - 1/6/02 - 4/10/11 (instrument died)
       % Note: filename dates are END of 3 day averaging period
       t1 = 37406;
       t2 = 40818;
       gg = [0 360 -90 90];
       dnm = [pth2 '3_day' slsh 'all' slsh];
     case 17
       % 7-day ave SSMI Microwave;  Global - 1/6/02 - 1/10/11 (instrument died)
       % Note: filename dates are END of 7 day averaging period
       t1 = 37406;
       t2 = 40815;
       gg = [0 360 -90 90];
       dnm = [pth2 'weekly' slsh];

     case {18,19,20,21,22}     
       % NOO-funded 2004 stitched archive. Australiasia, now kept updated
       % 1/10/93->
       itper = pr-17;
       t1 = t0 + [8674.5 8675.5 8677 8680 8681.5];
       t1 = t1(itper);
       %t2 = t0 + [12216.5 12215.5 12214 12212 12209.5];
       %t2 = t2(itper);
       t2 = inf;
       tper = [1 3 6 10 15];
       tper = tper(itper);       
       %gg = [79.987 190.03 -64.998 10.026];
       % Come in a fraction to avoid interp calc crashing right on boundaries
       gg = [79.99 190.02 -64.995 10.02];
       fnm = [pth3 'SSTcomp' num2str(tper) 'd_Aasia_'];
              
     case {23,24,25}     
       % Australiasia, NRT, daily files, lon .042, lat .036 deg, 0ct 2004 ->
       % NOTE: files are available prior to Oct 2004, at pth4b, but noone
       %       has asked for them yet.
       itper = pr-22;
       t1 = t0 + [12711 12710 12708];
       t1 = t1(itper);
       t2 = inf;
       tper = [3 6 10];
       tper = tper(itper);       
       %gg = [79.987 190.03 -64.998 10.026];
       gg = [79.99 190.02 -64.995 10.02];
       dnm = [pth4 'comp' num2str(tper) 'd' slsh 'Aasia' slsh];
     
     case {26,27}     
       % Australiasia, NRT, 0ct 2004 ->
       disp(['*** GET_SST_XY: Dataset ' num2str(pr) ' not yet generated']);       
       gg = [];       
     
     case 28
       % "Patchwork" V7. Australiasia, 1/6th degree lat and lon, 6 day space.
       t1 = 33768;
       t2 = 38214;
       gg = [84.916 185.083 -75.084 25.083];
       fnm = [pth5 'high_res_6_day_V7'];
       t0 = greg2time([1990 1 1 0 0 0]);
              
     case 29
       % "Patchwork" V8. Australiasia, 1/6th degree lat and lon, 6 day space.
       t1 = 38214;
       t2 = 38466;
       gg = [84.916 185.083 -75.084 25.083];
       fnm = [pth5 'high_res_6_day_V8'];
       t0 = greg2time([1990 1 1 0 0 0]);
              
     case 30
       %  Reynolds V7. Global, 1 degree lat and lon, 6 day space.
       t1 = 33768;
       t2 = 38214;
       gg = [0 360 -90 90];
       fnm = [pth5 'Reynolds_6_day_V7'];
       t0 = greg2time([1990 1 1 0 0 0]);
              
     case 31
       %  Reynolds V8 Global, 1 degree lat and lon, 6 day space.
       t1 = 38214;
       t2 = 38466;
       gg = [0 360 -90 90];
       fnm = [pth5 'Reynolds_6_day_V8'];
       t0 = greg2time([1990 1 1 0 0 0]);
              
     case 33
       %  Reynolds Weekly V2 Global, 1 degree lat and lon
       t1 = 29886;
       t2 = inf;
       gg = [0 360 -90 90];
       fnm = [pth6 'sst.wkmean_v2.'];
       t0 = greg2time([1800 1 1 0 0 0]);
              
     case 34
       %  NOAA Monthly V3 Global, 2 degree lat and lon
       % NOTE: the "time" variable contains month start dates, and
       % "time_bnds" values are month start and end, so suggests that values
       % actually relate to month average, ie centred on middle of month. I
       % have not yet allowed for this.    JRD 23/3/10       
       t1 = 19723;  % 1/1/1854
       t2 = inf;
       gg = [0 360 -88 88];
       fnm = [pth7 'sst.mnmean.'];
       t0 = greg2time([1800 1 1 0 0 0]);
              
     case 35
       %  Windsat daily V7 Global, 0.25 degree lat and lon
       %     JRD 21/2/2012       
       t1 = 40542;  % 1/1/2011  There is earlier data, but not yet downloaded
       t2 = inf;
       gg = [0 360 -89.9 89.9];
       greg = time2greg(tim);
       dnm = [pth8 num2str(greg(1)) '/'];
       t0 = greg2time([1981 1 1 0 0 0]);
              
     otherwise
       warning(['GET_SST_XY: do not know dataset ' num2str(pr)]);       
       gg = [];
   end

   jin = [];
   if ~isempty(gg)
      if tim>t1 & tim<t2 
	 xp = gg([1 1 2 2]); yp = gg([3 4 4 3]);
	 jin = jj(find(inpolygon(x(jj),y(jj),xp,yp)));
      end
   end

   vg = [];
   if isempty(jin)
      % no points to be found in this dataset
   elseif pr==1
      % PF OI 8km
      % Set up range of remaining required points, with border > grid
      % interval so that can interpolate to all points. Restrict to dataset
      % region to stop pfsst from complaining.
      rng = [min(x(jin))-.12 max(x(jin))+.12 min(y(jin))-.12 max(y(jin))+.12];
      rng([1 3]) = max([rng([1 3]); gg([1 3])]); 
      rng([2 4]) = min([rng([2 4]); gg([2 4])]); 
      [vg,glo,gla,timug] = pfsst(rng,time2greg(tim));
      timu(pr) = timug-julian([1900 1 1 0 0 0]);
   
   elseif any([3 5 6 7 10 11 12 13 14 15 16 17 23 24 25 35]==pr)
      % Daily composite netCDF files (previously accessed disimp versions)      
      % Get a directory listing; decode filenames to times; find nearest time
      dirl = dir([dnm '*.nc']);
      if ~isempty(dirl)	 
	 ftim = zeros([1 length(dirl)]);
	 idsh = findstr(dirl(1).name,'_');
	 if pr~=35 && length(idsh)~=ndsh(pr) 
	    disp([7 'Problem with filenames in ' dnm]);
	    dum = nan;
	 else
	    if isempty(idsh) || ndsh(pr)==0
	       idsh = 0;
	    else
	       idsh = idsh(ndsh(pr));
	    end
	    for ii = 1:length(dirl)
	       if dirl(ii).bytes > 0
		  yr = str2num(dirl(ii).name(idsh+(1:4)));
		  mo = str2num(dirl(ii).name(idsh+(5:6)));
		  da = str2num(dirl(ii).name(idsh+(7:8)));
		  ftim(ii) = greg2time([yr mo da 0 0 0]) - toff(pr);
	       else
		  ftim(ii) = nan;
	       end
	    end
	    [dum,itm] = min(abs(ftim-tim));
	 end
	 
	 if isnan(dum)
	    % no good files
	    vg = [];
	 else
	    fnm = [dnm dirl(itm).name];
	    if pr>12
	       tmp = getnc(fnm,'time');
	       if pr==35
		  % Average the 2 Windsat times and convert from secs to days
		  timu(pr) = (mean(tmp)/86400)+t0;
	       elseif ~isempty(tmp)
		  timu(pr) = tmp+t0 ;
	       end
	    else
	       timu(pr) = ftim(itm);	 
	    end
	    gla = getnc(fnm,'lat');
	    glo = getnc(fnm,'lon');
	    if any([16 17]==pr)
	       % Avoid missing values between 359.88E and 0.12E
	       glo(1) = -.01; glo(end) = 360.01;
	    end
	    if pr==35
	       % Windsat has orbits separated into ascending and descending
               % time arrays. I guess we use an average here where is
               % possible. Also convert from Kelvin to C.

	       switch passprf
		 case 0
		   i1=1; i2=2;
		 case 1
		   i1=2; i2=2;
		 case 2
		   i1=1; i2=1;
	       end
		
	       vg = getnc(fnm,'sea_surface_temperature',[i1 -1 -1],[i2 -1 -1]) ...
		    - 273.15;
	       flg = getnc(fnm,'confidence_flag',[i1 -1 -1],[i2 -1 -1]);
	       vg(flg~=0) = nan;
	       flg = getnc(fnm,'proximity_confidence',[i1 -1 -1],[i2 -1 -1]);
	       vg(flg<=2) = nan;

	       if passprf==0
		  vg2 = squeeze(vg(2,:,:));
		  vg =  squeeze(vg(1,:,:));
		  
		  jk = ~isnan(vg) & ~isnan(vg2);
		  vg(jk) = (vg(jk)+vg2(jk))./2;

		  jk = isnan(vg) & ~isnan(vg2);
		  vg(jk) = vg2(jk);
	       end
	    else
	       vg = getnc(fnm,'sst',[1 -1 -1],[1 -1 -1]);
	    end
	 end
      end
      
   %elseif pr==2 | pr==8 | pr==9
      % netCDF timeseries composite files
     % stim = getnc(fnm,'time')+t0;
     % [dum,itm] = min(abs(stim-tim));
     % timu(pr) = stim(itm);
     % gla = getnc(fnm,'lat');
     % glo = getnc(fnm,'lon');
     % vg = getnc(fnm,'sst',[itm -1 -1],[itm 1 1]);
   
   
   elseif any([18 19 20 21]==pr)
      % Yearly netCDF files
      fnm = sst_filnam(tim,fnm);
      stim = getnc(fnm,'time')+t0;
      [dum,itm] = min(abs(stim-tim));
      timu(pr) = stim(itm);      
      gla = getnc(fnm,'lat');
      glo = getnc(fnm,'lon');
      ix1 = sum(glo<min(x(jin)));
      ix2 = length(glo)+1 - sum(glo>max(x(jin)));
      iy1 = sum(gla>max(y(jin)));
      iy2 = length(gla)+1 - sum(gla<min(y(jin)));
      glo = glo(ix1:ix2);
      gla = gla(iy1:iy2);
      vg = getnc(fnm,'sst',[itm iy1 ix1],[itm iy2 ix2]);

      if ~isempty(mincl) && tim>greg2time([2005 1 1 0 0 0])
	 nclr = getnc(fnm,'nclear',[itm iy1 ix1],[itm iy2 ix2]);
      end
      
%      disp(fnm)
%      disp(num2str([itm iy1 ix1 iy2 ix2]))

   elseif pr==22
      % Yearly netCDF files with 15-day window, 6-day spacing, so interpolate
      % in time before interp in space.
      % Yearly netCDF files
      fnm1 = sst_filnam(tim,fnm);
      stim = getnc(fnm1,'time')+t0;
      nn = length(stim);

      gla = getnc(fnm1,'lat');
      glo = getnc(fnm1,'lon');      
      ix1 = sum(glo<min(x(jin)));
      ix2 = length(glo)+1 - sum(glo>max(x(jin)));
      iy1 = sum(gla>max(y(jin)));
      iy2 = length(gla)+1 - sum(gla<min(y(jin)));
      glo = glo(ix1:ix2);
      gla = gla(iy1:iy2);

      if tim<stim(1)
	 st2 = stim(1);
	 fnm2 = sst_filnam(tim,fnm,-1);	 
	 stim = getnc(fnm2,'time')+t0;
	 nn = length(stim);
	 st1 = stim(nn);
	 vg = getnc(fnm2,'sst',[nn iy1 ix1],[nn iy2 ix2]);
	 vg = cat(3,vg,getnc(fnm1,'sst',[1 iy1 ix1],[1 iy2 ix2]));
      elseif tim>stim(nn)
	 st1 = stim(nn);
	 fnm2 = sst_filnam(tim,fnm,1);	 
	 stim = getnc(fnm2,'time')+t0;
	 st2 = stim(1);
	 vg = getnc(fnm1,'sst',[nn iy1 ix1],[nn iy2 ix2]);
	 vg = cat(3,vg,getnc(fnm2,'sst',[1 iy1 ix1],[1 iy2 ix2]));
      else
	 it1 = find(stim<tim);
	 it1 = it1(end);
	 st1 = stim(it1);
	 st2 = stim(it1+1);	 
	 vg = getnc(fnm1,'sst',[it1 iy1 ix1],[it1+1 iy2 ix2]);
	 vg = shiftdim(vg,1);
      end
      timu(pr) = tim;  

      trat = (tim-st1)/(st2-st1);
      vg = squeeze(vg(:,:,1)) + (trat.*squeeze(vg(:,:,2)-vg(:,:,1)));
   
   elseif any([28 29 30 31]==pr)
      % Single netCDF file with 6-day spacing, so interpolate in time
      stim = getnc(fnm,'time')+t0;
      
      if isempty(stim)
	 disp(['GET_SST_XY, dataset ' num2str(pr) ' This file might be gzipped?']);
      else
	 gla = getnc(fnm,'latitude');
	 glo = getnc(fnm,'longitude');      
	 if pr==28 || pr==29
	    iy1 = sum(gla<min(y(jin)));
	    iy2 = length(gla)+1 - sum(gla>max(y(jin)));
	 else
	    % Reynolds files have descending latitude and lon = [.5 to 359.5]
	    glo(1) = -.01; glo(end) = 360.01;
	    iy1 = sum(gla>max(y(jin)));
	    iy2 = length(gla)+1 - sum(gla<min(y(jin)));
	 end
	 ix1 = sum(glo<min(x(jin)));
	 ix2 = length(glo)+1 - sum(glo>max(x(jin)));
	 glo = glo(ix1:ix2);
	 gla = gla(iy1:iy2);

	 it1 = find(stim<tim);
	 it1 = it1(end);
	 vg = getnc(fnm,'sst',[it1 iy1 ix1],[it1+1 iy2 ix2]);

	 stm = stim([it1 it1+1]);
	 trat = (tim-stm(1))/diff(stm);
	 vg = sq(vg(1,:,:)) + (trat.*sq(vg(2,:,:)-vg(1,:,:)));
      
	 timu(pr) = tim;  
      end
      
   elseif pr==33
      % Two netCDF file with 7-day spacing, so interpolate in time
      % time=days since 1800. Note that variable 'time' is the start of
      % weekly averaging period - so we use mid-point of "time_bnds" instead
      % to give us an 'estimate centre time'.
      if tim<greg2time([1990 1 1 0 0 0])
	 fnm = [fnm '1981-1989'];
      else
	 fnm = [fnm '1990-present'];
      end
      stim = getnc(fnm,'time_bnds');
      stim = t0 + (stim(:,1)+stim(:,2))./2;

      gla = getnc(fnm,'lat');
      glo = getnc(fnm,'lon');      
      % A cludge to retrieve values for locations between 359.5E and 0.5E
      glo(1) = -.01; glo(end) = 360.01;
      ix1 = sum(glo<min(x(jin)));
      ix2 = length(glo)+1 - sum(glo>max(x(jin)));
      % Reynolds folk have latitude +89.5 down to -89.5
      iy1 = sum(gla>max(y(jin)));
      iy2 = length(gla)+1 - sum(gla<min(y(jin)));
      glo = glo(ix1:ix2);
      gla = gla(iy1:iy2);

      it1 = find(stim<tim);
      it1 = it1(end);

      % Reynolds file has incorrect 'valid_range' attribute which causes
      % most data to be trashed, so have to do our own scaling!
      ncid = netcdf([fnm '.nc'],'nowrite');
      scf = ncid{'sst'}.scale_factor(:);
      ado = ncid{'sst'}.add_offset(:);
      msv = ncid{'sst'}.missing_value(:);
      vg = ncid{'sst'}(it1:it1+1,iy1:iy2,ix1:ix2);
      ncid = close(ncid);

      ims = find(vg==msv);
      vg(ims) = nan;
      vg = (vg*scf)+ado;

      %vg = getnc(fnm,'sst',[it1 iy1 ix1],[it1+1 iy2 ix2]);

      stm = stim([it1 it1+1]);
      trat = (tim-stm(1))/diff(stm);
      vg = sq(vg(1,:,:)) + (trat.*sq(vg(2,:,:)-vg(1,:,:)));
      
      timu(pr) = tim;  
   end
   
   % Have some gridded data, so interpolate to our locations and set dset 
   % appropriately where data is obtained at our locations.
   if ~isempty(vg)
      if ~isempty(nclr)
	 vg(nclr<mincl) = NaN;
	 nclr = [];
      end
      if pr==35
	 % Windsat longitude is [-180 180] 
	 jk = x>180;
	 x(jk) = x(jk)-360;
      end
      if sprd
	 vg = spread1cell(vg);
      end
      sst(jin) = interp2(glo,gla,vg,x(jin),y(jin));
      dset(jin) = ~isnan(sst(jin)).*pr;      
   end
   
   % Find what locations still need filling, and test whether to continue
   jj = find(isnan(sst));
   req = ~isempty(jj) & ip<length(pref) & ~(sngl & any(~isnan(sst(:))));
   ip = ip+1;
end

% If require only first timestamp, then organise it.
if tmu1
   ii = find(timu);
   if isempty(ii)
      timu = 0;
   else
      timu = timu(ii(1));
   end
end

%---------------------------------------------------------------------------
