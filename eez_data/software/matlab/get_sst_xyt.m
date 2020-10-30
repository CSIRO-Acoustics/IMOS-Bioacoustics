% GET_SST_XYT  Retrieve data to given times and locations from the 2004
%   4km Stitched SST Archive.   Note: only actually interpolates for 15-day 
%   composites as these are at 6-day spacing. 10-day composites are at 2-day
%   spacing, and all the rest at 1-day, so we just take nearest time slice
%   for all of these.
%
% NOTE:  Uses Matlab interp functions which give a Nan if any adjacent 
%        input value is Nan, so data recovery rates can be poor for
%        small-time-window composites (which have more gaps.)
%
%  See also GET_SST, which returns fuul-res SST dataset grid in given region,
%  and GET_SST_XY, which interpolates in space but evaluates all points at 
%  just one time value.
%
% INPUT
%  x,y  - vector or matrices of locations
%  tim  - Unix time for each location (can use greg2time)
%  itper- single code for composite time-window size:
%        1=1   2=3  3=6day  4=10day  5=15day
%
% OUTPUT
%  sst  - SST at x,y locations (nan where no data available)
%
% SEE ALSO: ~ridgway/matlab/altim/get_sstPW_xyt.m  for access to "Patchwork"
%
% SEE http://www.marine.csiro.au/remotesensing/oceancurrents/  OR
%     http://www.marine.csiro.au/eez_data/doc/sst_datasets.html
%
% Jeff Dunn  CSIRO CMAR 10/6/04, 16/8/07
%
% USAGE: sst = get_sst_xyt(x,y,tim,itper)

function sst = get_sst_xyt(x,y,tim,itper)

% ABOUT Algorithm
%
% Have tried many ways with this code, and each is very slow for some given
% situation. For sparse data extracting a chunk spanning the whole time of
% a file, and interpolating in that, is much more costly than doing lots
% of little extractions for just the times required. This present coding is
% a compromise, with some attempt to keep the code simple but fast.

sst = repmat(nan,size(x));

t70 = greg2time([1970 1 1 0 0 0]);
tperiod = [1 3 6 10 15];
tper = tperiod(itper);
tdel = [1 1 1 2 6];
tdel = tdel(itper);
tstart = [34241.5 34242.5 34244 34247 34248.5];
tstart = tstart(itper);

fnm = path_pc_or_nix('imgjj/sensor/avhrr/sstcr04/yearfiles/');
fnm = [fnm 'SSTcomp' num2str(tper) 'd_Aasia_'];

la = getnc([fnm '1994'],'lat');
lo = getnc([fnm '1994'],'lon');

kk = find(x>min(lo) & x<max(lo) & y>min(la) & y<max(la) & tim>tstart);

% x spacing .042 degree, y .036 degree, y stored in descending order!
nlo = length(lo);
ix = 1+(x-min(lo))./.042;
nla = length(la);
iy = nla+1-((y-min(la))./.036);


fnam = sst_filnam(min(tim(kk)),fnm);

while ~isempty(kk) && exist([fnam '.nc'],'file')

   ftim = getnc(fnam,'time') + t70;
   ntim = length(ftim);
   
   itm = 1+((tim(kk)-ftim(1))./tdel);
   
   if tdel > 2
      % 15 day window composites at 6 day spacings: interpolate time
      
      jj = find(itm<ntim);
      itm = itm(jj);
      kj = kk(jj);
   

      for ll = unique(floor(itm(:)'))
	 mm = find(itm>=ll & itm<ll+1);
	 
	 xs = ix(kj(mm));
	 ys = iy(kj(mm));	    
	 x0 = max([1 floor(min(xs))]);
	 y0 = max([1 floor(min(ys))]);
	 x1 = min([nlo ceil(max(xs))]);
	 y1 = min([nla ceil(max(ys))]);
	 
	 if ll == 0
	    % Times in gap between this and preceding file
	    fnaml = sst_filnam(min(tim(kj))-2*tdel,fnm);
	    ltim = getnc(fnaml,'time') + t70;
	    nltm = length(ltim);
	    sstl = getnc(fnaml,'sst',[nltm y0 x0],[nltm y1 x1]);
	    sstf = getnc(fnam,'sst',[1 y0 x0],[1 y1 x1]);
	    sstf = cat(1,shiftdim(sstl,-1),shiftdim(sstf,-1));
	 else
	    sstf = getnc(fnam,'sst',[ll y0 x0],[ll+1 y1 x1]);
	 end
	 
	 sst(kj(mm)) = interp3(sstf,1+ys-y0,1+rem(itm(mm),1),1+xs-x0);
      end
   
   else
      % 1 or 2 day spaced estimates - just select nearest time slice.
	 
      itm = round(itm);
      jj = find(itm>=1 & itm<=ntim);
      itm = itm(jj);
      kj = kk(jj);
	 
      for ll = unique(itm(:)')
	 mm = find(itm==ll);
	 if ~isempty(mm)
	    xs = ix(kj(mm));
	    ys = iy(kj(mm));	    
	    x0 = max([1 floor(min(xs))]);
	    y0 = max([1 floor(min(ys))]);
	    x1 = min([nlo ceil(max(xs))]);
	    y1 = min([nla ceil(max(ys))]);
	    sstf = getnc(fnam,'sst',[ll y0 x0],[ll y1 x1]);
	    sst(kj(mm)) = interp2(sstf,1+xs-x0,1+ys-y0);
	 end
      end   
   end

   kk(jj) = [];
   if ~isempty(kk)
      fnam = sst_filnam(min(tim(kk))+2*tdel,fnm);
   end
end      

%---------------------------------------------------------------------------
