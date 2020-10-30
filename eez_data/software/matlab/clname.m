% CLNAME: Get name of climatology file
% INPUT: 
%    property - property name (eg 'salinity' or just 's')
%    pth    - path to netCDF map file
%    fname   - input file name component, if non-CARS. Full name is built as
%              follows:    filename = [fullproperty '_' fname '.nc']
%            eg: property='s', fname='maps'  then filename='salinity_maps.nc'
%
% OUTPUT: 
%    clnam   Constructed climatology name (not including '.nc')
%    ncf     Optional: open the file; returns the netcdf toolbox object
%
% MODS: 6/8/10 Now default to CARS_latest
%
% USAGE: [clnam,ncf] = clname(property,pth,fname)

function [clnam,ncf] = clname(property,pth,fname)

persistent CLNAME_def_warn

if length(property)==1
  property = [property ' '];
end

if strcmp(property(1),'t')
  property = 'temperature';
elseif strcmp(property(1),'o')
  property = 'oxygen';
elseif strcmp(property(1),'n')
  property = 'nitrate';
elseif strcmp(property(1),'p')
  property = 'phosphate';
elseif strcmp(property(1:2),'si')
  property = 'silicate';
elseif strcmp(property(1),'s')
  property = 'salinity';
end
  
if nargin<2 | isempty(pth)
   pth = path_pc_or_nix('climatologies/CARS/');
end

if nargin<3 | isempty(fname)
   if isempty(CLNAME_def_warn)
      CLNAME_def_warn = 1;
      disp('CLNAME - "CARS_latest" version being used by default.')
      disp('[  This warning once per session only  ]')
   end
  fname='CARS_latest';
end

clnam = [pth property '_' fname];
if ~exist([clnam '.nc'],'file')
   pth = path_pc_or_nix('eez_data/atlas/');
   clnam2 = [pth property '_' fname];
   if ~exist([clnam2 '.nc'],'file')
      disp(['*** Cannot find a CARS file called ' clnam '.nc']);
      disp('Using latest CARS version instead');
      fname='CARS_latest';
      clnam = [pth property '_' fname];
   else
      clnam = clnam2;
   end
end
   
if nargout>1
  ncf = netcdf([clnam '.nc'],'nowrite');
  if isempty(ncf)
    disp(['CLNAME: *** Cannot open CARS file ' clnam]);
  end
end

% -------------------------------------------------------
