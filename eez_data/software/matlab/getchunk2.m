% GETCHUNK2: Extract a 3D chunk of a netCDF file map.
%            Same functionality as getchunk but coded to suit new Matlab
%
% INPUT:
%    property - property name (eg 'salinity' or just 's')
%  Optional....
%    deps    - depths (in +ve m, sorted, no duplicates) of levels to extract
%    region  - [w e s n] boundaries of maps to extract
%    pth     - path to netCDF map file (including trailing slash)              
%    fname   - input file name component, if not 'cars2000'. Full name is built as
%              follows:    filename = [fullproperty '_' fname '.nc']
%                eg: property='s', fname='v1'  then filename='salinity_v1.nc'
%    opt       1= [default] return lat/lon for x,y
%              2= return grid indices (ix,iy) for x,y
%              -1 or -2: as above, but return empty args for all except mn,x,y
%    vart      1=seasonal or mean  9=SD  10=RMSMR  11=RMSR  [def 1]
%
% OUTPUT:
%    mn     - [nlat,nlon,ndep] mean (or vart-specified variable) field
% Optional....
%    an     - [nlat,nlon,ndep] complex annual harmonics, [] if none available
%    sa     - [nlat,nlon,ndep] complex semiannual harmonics, "  "  "   " 
%    rq     - [nlat,nlon,ndep] data source radius (km)
%    dets   - [ndep,det_len] mapping details text strings
%                opt = 1                         opt = 2    [see 'opt' above])
%    x      - [nlat,nlon] longitudes         [nlon] map grid coords
%    y      - [nlat,nlon] latitudes          [nlat] map grid coords
%
% * If ONLY ONE DEPTH specified or available, 3D outputs collapse to 2D.
%
% Copyright (C) J R Dunn, CSIRO Marine Research  Aug 2012
%
% Eg:  [mn,a1,a2,a3,a4,x,y]= getchunk2('t',[0 10 20 50],[90 100 -40 0],[],[],1);
% 
% USAGE: [mn,an,sa,rq,dets,x,y] = getchunk2(property,deps,region,pth,fname,opt,vart)

function [mn,an,sa,rq,dets,x,y] = getchunk2(property,deps,region,pth,fname,opt,vart)

[mn,an,sa,rq,dets,x,y] = deal([]);

if nargin<7 || isempty(vart)
   vart = 1;
end
if nargin<6 || isempty(opt)
  opt = 1;
end
if nargin<5 
   fname = [];
end
if nargin<4
   pth = [];
end

cfnm = clname(property,pth,fname);
cfnm = [cfnm '.nc'];

if nargin<3; region=[]; end

if vart>1
   % If asking for variability fields, which never have corresponding
   % seasonal fields, then silly to return seasonal value fields.
   opt = -(abs(opt));
end

getan = 0; an = [];
getsa = 0; sa = [];

% Check existence of any requested data

if ~exist(cfnm,'file')
  error(['Cant open file ' cfnm]);
end
   
% Must allow for non-contiguous vector of depths

depths = round(ncread(cfnm,'depth'));
if nargin<2 | isempty(deps)
   deps = depths;
else
   orig_dep = deps;
   deps = round(deps(:)');
end
ndep = length(deps);

[tmp,idin,idout] = intersect(depths,deps);
if length(idout)<ndep
   jj = 1:ndep;
   jj(idout) = [];
   if isempty(idout)
      disp('None of these depths present in this file - returning');
      return
   else
      disp(['These depths are not present in this file: ' num2str(orig_dep(jj))]);
   end
end

dep1 = idin(1);
idin = 1+idin-dep1;
depn = idin(end);

info = ncinfo(cfnm);
dimnm{length(info.Dimensions)} = [];
[dimnm{:}] = deal(info.Dimensions.Name);
varnm{length(info.Variables)} = [];
[varnm{:}] = deal(info.Variables.Name);

if nargout>1 & opt>=0
  % Find out at what depths these variables exist  
  iid = strmatch('depth_ann',dimnm);
  pre06 = isempty(iid);
  if pre06
     % pre 2006 versions
     iid = strmatch('depth_timefit',dimnm);
  end
     
  if isempty(iid)
     disp('No temporal info.');
  else     
     dtmax = info.Dimensions(iid).Length;
     [tmp,itin,ito1] = intersect(depths(1:dtmax),deps);
     if ~isempty(itin)
	tep1 = itin(1);
	tepn = itin(end);
	itin = 1+itin-(tep1);
	getan = 1;
	iiv = strmatch('sa_cos',varnm);
	if nargout>2 && ~isempty(iiv)
	   if pre06
	      getsa = 1;
	      tep2 = tep1;
	      tepn2 = tepn;
	      ito2 = ito1;
	      itin2 = itin;	      
	   else
	      iid = strmatch('depth_semiann',dimnm);
	      dmax2 = info.Dimensions(iid).Length;
	      [tmp,itin2,ito2] = intersect(depths(1:dmax2),deps);
	      if ~isempty(itin2)
		 tep2 = itin2(1);
		 tepn2 = itin2(end);
		 itin2 = 1+itin2-(tep2);
		 getsa = 1;
	      end
	   end
	end
     end
  end
else
  rq = [];
  dets = [];
end


iiv = strmatch('map_details',varnm);
got_dstr = ~isempty(iiv);

iiv = strmatch('radius_q',varnm);
got_rq = ~isempty(iiv);


Xg = ncread(cfnm,'lon');
Yg = ncread(cfnm,'lat');
if min(size(Xg))==1
   [Xg,Yg] = meshgrid(Xg,Yg);
end


if isempty(region)
   x1 = 1; y1 = 1;
   x2 = length(Xg(1,:));
   y2 = length(Yg(:,1));
else
   % Note that even for rotated grids this will find the minimum 
   % grid-rectangular area enclosing all of 'region'
   [iy,ix] = find(Yg>=region(3)& Yg<=region(4)& Xg>=region(1)& Xg<=region(2));
   y1 = min(iy); y2 = max(iy);
   x1 = min(ix); x2 = max(ix);
end

ny = 1+y2-y1;
nx = 1+x2-x1;
  

% getnc outputs are 3D unless the depth dimension is scalar.
% Get mean and shape it. DIMS is a function at the end of this file

varnms = {'mean','','','','','','','','RMSspatialresid','RMSresid','std_dev'};
varn = varnms{vart};

tmp = ncread(cfnm,varn,[x1 y1 dep1],[nx ny depn]);
if depn==1
   mn = tmp';
else
   mn = nan([ny nx ndep]);
   for ii = 1:ndep
      mn(:,:,idout(ii)) = squeeze(tmp(:,:,idin(ii)))';
   end
end
mn = squeeze(mn);


if getan
   vcos = ncread(cfnm,'an_cos',[x1 y1 tep1],[nx ny 1+tepn-tep1]);
   vsin = ncread(cfnm,'an_sin',[x1 y1 tep1],[nx ny 1+tepn-tep1]);
   if ndims(vcos)==2
      an = vcos' + (vsin').*1i;
   else
      an = nan([ny nx ndep]);
      for ii = 1:length(itin)
	 an(:,:,idout(ii)) = vcos(:,:,itin(ii))' + vsin(:,:,itin(ii))'.*1i;
      end  
   end
end

if getsa
   vcos = ncread(cfnm,'sa_cos',[x1 y1 tep2],[nx ny 1+tepn2-tep2]);
   vsin = ncread(cfnm,'sa_sin',[x1 y1 tep2],[nx ny 1+tepn2-tep2]);
   if ndims(vcos)==2
      sa = vcos' + (vsin').*1i;
   else
      sa = nan([ny nx ndep]);
      for ii = 1:length(itin2)
	 sa(:,:,idout(ii)) = vcos(:,:,itin2(ii))' + vsin(:,:,itin2(ii))'.*1i;
      end  
   end
end

if nargout>3 && opt>=0
   if ~got_rq
      rq = zeros(size(mn));
  else
     tmp = ncread(cfnm,'radius_q',[x1 y1 dep1],[nx ny depn]);
     if ndims(tmp)==2
	rq = tmp';
     else
	rq = nan([ny nx ndep]);
	for ii = 1:length(idin)
	   rq(:,:,idout(ii)) = squeeze(tmp(:,:,idin(ii)))';
	end
     end
  end
end
  

if nargout>4 && opt>=0
   if got_dstr
      dets = ncread(cfnm,'map_details',[dep1 1],[depn inf]);
      dets(idout,:) = dets(idin,:);
   else
      dets = repmat(' ',[ndep 1]);
   end
end

if nargout > 5
   if abs(opt)==2
      y = y1:y2;
      x = x1:x2;
   else
      x = Xg(y1:y2,x1:x2);
      y = Yg(y1:y2,x1:x2);
   end
end

% --------------------------------------------------------------
