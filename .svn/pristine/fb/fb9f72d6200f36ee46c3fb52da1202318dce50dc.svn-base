% GET_BATH  Get ocean depth at given locations, using a range of datasets.
%
%  **** See www.marine.csiro.au/eez_data/doc/bathy.html 
%
%  **** NOTE Default dset changed from [5 6] to [5 10] on 6/5/09
%
% INPUTS
%  lon,lat  Locations at which depth is required
%  dset   vector of datasets in preferred order    Default:  [5 10]
%         1  Terrainbase  
%         2  NGDCv8.2  
%         3  AusBath15  [redundant]   
%         4  AGSO98  
%         5  AGSO2002
%         6  GEBCO'03   [superceded by 10]
%         7  ETOPO2v2   
%            (NOTE:  8 - "GA 9sec" is accessed only via GA9SECBATH) 
%         9  S2004
%        10  GEBCO 1min v2 (2008)
%        11  GEBCO 30sec (2008)
%  rpt    [redundant]
%
% OUTPUTS
%  deps   depths (m), -ve downwards, NaN where no value, interpolated to
%         locations EXCEPT for datasets 1,2,7, which are just sampled from 
%         datasets.
%
% SUPERCEDES  get_bath15, agso_bath_xy
%
% Jeff Dunn   CSIRO Marine Research   8/1/2003 - 27/7/07
%
% SEE ALSO   get_bath_agso.m   (to get AGSO 2002 at full resolution)
%
% USAGE: deps = get_bath(lon,lat,dset);

function deps = get_bath(x,y,dset,dummy)

% Mods: 9/5/03 Pre-check for presence of datasets, and resort to others
%       if some are missing.
%       27/7/07 Add dataset 7, remove redundant "rpt" argument from comments
%       26/11/07  Add dataset 9
%       23/4/09  Add datset 10,11
%       6/5/09  Change default to [5 10]. Remove last trace of "rpt"
%       26/7/10  Chunksize increase from 500->5000. (idea from Paul Durack)
%       5/4/12 Mod to allow interpolation to values at 0 or 360E
%       2/2/13 Chunksize dropped back to 1000, although crash was prob due to
%              error in the 26/7/10 change.
%	21/8/14 Changed path for /home/eez_data as it's been moved. (RAS)

ncquiet;
   
if nargin<3 | isempty(dset)
   dset = [5 10];
elseif any(dset==8)
   disp('Use ga9secbath.m to access dataset 8')
   dset(dset==8) = [];
   if isempty(dset)
      deps = [];
      return
   end
end

igd = 0;
dtmp = [];
fnms{1} = [];
for ii = dset(:)'
   switch ii
     case 1
       fnm = path_pc_or_nix('eez_data/bath/terrainbase.mat');
       aa = exist(fnm,'file');
     case 2
       fnm = path_pc_or_nix('netcdf-data/topo_ngdc_8.2.nc');
       aa = exist(fnm,'file');
     case 3
       fnm = path_pc_or_nix('eez_data/bath/ausbath15.mat');
       aa = exist(fnm,'file');
       if ~aa
	  fnm = path_pc_or_nix('dunn/bath/ausbath15.mat');
	  aa = exist(fnm,'file');
       end
     case 4
       fnm = path_pc_or_nix('netcdf-data/bath_agso_98.nc');
       aa = exist(fnm,'file');       
     case 5
       fnm = path_pc_or_nix('netcdf-data/bath_agso_2002.nc');
       aa = exist(fnm,'file');
     case 6
       fnm = path_pc_or_nix('netcdf-data/gebco_2003.nc');
       aa = exist(fnm,'file');
     case 7
       fnm = path_pc_or_nix('netcdf-data/ETOPO2v2c_f4_cmar.nc');
       aa = exist(fnm,'file');
     case 9
       fnm = path_pc_or_nix('netcdf-data/bath_s2004.nc');
       aa = exist(fnm,'file');
     case 10
       fnm = path_pc_or_nix('eez_data/bath/gebco08_1min.nc');
       aa = exist(fnm,'file');
     case 11
       fnm = path_pc_or_nix('eez_data/bath/gebco08_30sec.nc');
       aa = exist(fnm,'file');
     otherwise
       aa = 0;
   end
   if aa>0
      igd = igd+1;
      dtmp(igd) = ii;
      fnms{igd} = fnm;
   else
      disp(['GET_BATH: dataset ' num2str(ii) ' not used because gone missing.'])
   end
end
if igd==0
   disp('GET_BATH: Cannot find any of the bathy files - no depths returned.')
   return
end
dset = dtmp;


deps = repmat(nan,size(x));

ii = 1:prod(size(x));

ids = 1;
while ~isempty(ii) & ids<=igd
   switch dset(ids)
     case 1
       % Terrainbase
       deps(ii) = topo(y(ii),x(ii),fnms{ids});
       ii = [];
     
     case 2
       % NGDC
       deps(ii) = topongdc(y(ii),x(ii),fnms{ids});
       ii = ii(find(y(ii)<-72 | y(ii)>72));
     
     case 3
       % AusBath15
       jj = find(x(ii)>=109 & x(ii)<=156 & y(ii)>=-45 & y(ii)<=-1);

       if ~isempty(jj)
	  ji = ii(jj);
	  if ~exist('AusBath15','var')
	     persistent AusBath15
	     load(fnms{ids});
	  end
	  
	  xx = 1+((x(ji)-109)*15);
	  yy = 1+((-1-y(ji))*15);
	  deps(ji) = interp2(AusBath15,xx,yy,'*linear');
	  ii(jj) = [];
       end

     case {4,5,6,9,10,11}
       % AGSO 98, AGSO 2002, GEBCO 2003 1min, GEBCO 1min v2, GEBCO08 30sec
       if dset(ids)==10
	  lo = 0:(1/60):360;
	  la = -90:(1/60):90;
       elseif dset(ids)==9 | dset(ids)==11
	  lo = getnc(fnms{ids},'longitude');
	  la = getnc(fnms{ids},'latitude');
       else
	  lo = getnc(fnms{ids},'lon');
	  la = getnc(fnms{ids},'lat');
       end

       jj = find(x(ii)>=min(lo) & x(ii)<=max(lo) & y(ii)>=min(la) & y(ii)<=max(la));
       ji = ii(jj);

       if ~isempty(ji)
	  % Broaden region slightly (by .1 degree) so extracted chunk encloses
	  % all points
	  ix = find(lo>=min(x(ji))-.1 & lo<=max(x(ji))+.1);
	  iy = find(la>=min(y(ji))-.1 & la<=max(y(ji))+.1);

	  lo(lo==0) = -.02;    % Enable interpolation to x=0
	  lo(lo==360) = 360.02;    % Enable interpolation to x=360

	  % If a large region required, break it into chunks to avoid causing
	  % a crash due to memory overload. 
	  if length(iy)*length(ix) > 1000000
	     %chsz = 5000;    % Updated from 500 on 26/7/10
	     chsz = 1000;    % Changed back to 1000 / 1000000 on 2/12/13
                             % because too much for some machines
	     ixch = min(ix):chsz:max(ix)+chsz;
	     if ixch(end)>length(lo); ixch(end)=length(lo); end
	     iych = min(iy):chsz:max(iy)+chsz;
	     if iych(end)>length(la); iych(end)=length(la); end

	     for kx = 1:(length(ixch)-1)
		loc = lo(ixch(kx):ixch(kx+1));
		for ky = 1:(length(iych)-1)
		   lac = la(iych(ky):iych(ky+1));
		   ki = find(x(ji)>=min(loc) & x(ji)<max(loc) & y(ji)>=min(lac) & y(ji)<max(lac));
		   if ~isempty(ki)
		      kj = ji(ki);
		      ji(ki) = [];
		      dg = getnc(fnms{ids},'height',[iych(ky) ixch(kx)],[iych(ky+1) ixch(kx+1)]);
		      if length(loc)==1	  
			 % Degenerate case where only want points on boundary of dataset
			 deps(kj) = interp1(lac,dg,y(kj));
		      elseif length(lac)==1	  
			 % Ditto
			 deps(kj) = interp1(loc,dg,x(kj));
		      else
			 [xc,yc] = meshgrid(loc,lac);
			 deps(kj) = interp2(xc,yc,dg,x(kj),y(kj));
		      end
		   end
		end
	     end
	     
	  else
	     % Small region - don't need to break into chunks
	     dg = getnc(fnms{ids},'height',[iy(1) ix(1)],[iy(end) ix(end)]);
	     [lo,la] = meshgrid(lo(ix),la(iy));
	     if length(ix)==1	  
		% Degenerate case where only want points on boundary of dataset
		deps(ji) = interp1(la,dg,y(ji));
	     elseif length(iy)==1	  
		% Ditto
		deps(ji) = interp1(lo,dg,x(ji));
	     else
		deps(ji) = interp2(lo,la,dg,x(ji),y(ji));
	     end
	     
	  end

	  % Remove from the list only points for which we have obtained data.
	  ll = find(~isnan(deps(ii(jj))));
	  ii(jj(ll)) = [];
	  
	  clear jj ji lo la dg xc yc
       end
       
     case 7
       % ETOPO2v2
       deps(ii) = etopo2v2(y(ii),x(ii));
       ii = [];
     
     otherwise
       disp(['GET_BATH: Do not understand dataset ' num2str(dset(ids))]);
   end   
   ids = ids+1;
end
   
%---------------------------------------------------------------------------
