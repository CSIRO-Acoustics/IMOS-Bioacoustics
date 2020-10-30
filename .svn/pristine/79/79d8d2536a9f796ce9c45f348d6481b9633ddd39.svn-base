% GET_WILLIS  The Willis QC-ed T dataset is split into regional files so
%   you do not have to load the whole thing if you just want a small region.
%   This function just concatonates data from the required regional files.
%
% INPUT  src - can be empty (all data), or [w e s n], or polygon definition,
%              or one of the 5 subset regions by name  [default: all]
%        deps - depth levels (in metres) to extract  [default: all depths]
%        tao  - 0=no TAO  1=TAO only  2=TAO plus all other data  [0]
%
% OUTPUT 
%
% MODS:  6/11/07 Separate TAO into another set of files, and set up access
%                selection.
%
% USAGE: [lo,la,tim,tz,botdep,ptyp,pindx,styp,sindx,cpn] = get_willis(src,deps,tao)

function [lo,la,tim,tz,botdep,ptyp,pindx,styp,sindx,cpn] = get_willis(src,deps,tao)

if nargin==0
   src = [];
end
if nargin<3 || isempty(tao)
   tao = 0;
end

wdep = 0:10:750;
if nargin<2 | isempty(deps)
   ideps = 1:length(wdep);
else
   ideps = find(ismember(wdep,deps));
   if length(ideps)<length(deps)
      error('Willis has only depths 0:10:750')   
   end
end

lo = []; la = []; tim = []; tz = []; botdep = []; cpn = [];
ptyp = []; pindx = []; styp = []; sindx = []; cflg = [];

pth = path_pc_or_nix('eez_data/willis/');
if tao==0
   fnm = {'north_2004','sth_atl_2004','west_2004','central_2004','east_2004'};	
   regset = {1,2,3,4,5};
elseif tao==1
   fnm = {'west_2004_tao','central_2004_tao','east_2004_tao'};	
   regset = {[],[],3,4,5};
elseif tao==2
   fnm = {'north_2004','sth_atl_2004','west_2004','central_2004','east_2004',...
	 'west_2004_tao','central_2004_tao','east_2004_tao'};
   regset = {1,2,[3 6],[4 7],[5 8]};
end
nset = length(fnm);

jset = [];
if isempty(src)
   jset = 1:nset;
elseif max(size(src))==1
   if ~ismember(src,1:nset)
      error(['GET_WILLIS: src unknown: ' num2str(src)])
   else
      jset = src;
   end
elseif min(size(src))==1
   if src(4)>30; jset = regset{1}; end
   if src(3)<30
      if src(1)<20 | src(2)>290; jset = [jset regset{2}]; end
      if src(1)<140 & src(2)>20; jset = [jset regset{3}]; end
      if src(1)<210 & src(2)>140; jset = [jset regset{4}]; end
      if src(1)<290 & src(2)>210; jset = [jset regset{5}]; end
   end
else
   if any(src(:,2)>30); jset = regset{1}; end
   if any(src(:,2)<30)
      % This is highly imperfect, but it is tricky!
      if any(src(:,1)<20 | src(:,1)>290); jset = [jset regset{2}]; end
      if any(src(:,1)<140 & src(:,1)>20); jset = [jset regset{3}]; end
      if any(src(:,1)<210 & src(:,1)>140); jset = [jset regset{4}]; end
      if any(src(:,1)<290 & src(:,1)>210); jset = [jset regset{5}]; end
   end
end

for ii = jset
   load([pth fnm{ii}]);
   lo = [lo; lon];
   la = [la; lat];
   tim = [tim; time];
   tz = [tz; t(:,ideps)];
   botdep = [botdep; bdep];
   ptyp = ptypes;
   pindx = [pindx; ptindx];
   styp = source;
   sindx = [sindx; src_indx];
   cpn = [cpn; stnno];
   cflg = [cflg; t_castflag];
end

if isempty(src) | max(size(src))==1
   % return everything
else
   if min(size(src))==1
      ii = find(cflg==0 & lo>=src(1) & lo<=src(2) & la>=src(3) & la<=src(4));
   else
      ii = find(cflg==0 & inpolygon(lo,la,src(:,1),src(:,2))); 
   end
   lo = lo(ii);
   la = la(ii);
   tim = tim(ii);
   tz = tz(ii,:);
   botdep = botdep(ii);
   pindx = pindx(ii);
   sindx = sindx(ii);
   cpn = cpn(ii);   
end

%--------------------------------------------------------------------------
