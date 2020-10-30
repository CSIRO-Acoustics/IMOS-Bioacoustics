% Get the "cars_version" global attribute from a netCDF file
%
% INPUT  
%    fnm   full path name [if non-CARS netCDF file]
%    cnm   CARS name part [if CARS netCDF file]
%    pth   path [if CARS file and non-standard path]
%
% eg  cver = cars_vers_name('myfile.nc');
% eg  cver = cars_vers_name([],'cars2005a');
%
% NOTE: cver=[] if global attribute 'cars_version' is missing from the file
%
% USAGE: cver = cars_vers_name(fnm,cnm,pth);

function cver = cars_vers_name(fnm,cnm,pth)

cver = [];

if nargin<3
   pth = [];
end

if isempty(fnm)
   % Try different properties until we find a file, in case this version is 
   % not available for all properties.
   pnm = {'t','s','o','si','p','n'};
   ii = 1;
   while ii<=6      
      fnm = clname(pnm{ii},pth,cnm);
      fnm = [fnm '.nc'];
      if exist(fnm,'file')
	 ii = 99;
      elseif ii==6
	 % Didn't find a file
	 disp(['No files available with name ' cnm])
	 return
      else
	 ii = ii+1;
      end
   end
else
   nn = length(fnm);
   if ~strcmp(fnm((nn-2):nn),'.nc')
      fnm = [fnm '.nc'];
   end
end

[attval,attname] = attnc(fnm);

ii = strmatch('cars_version',attname,'exact');
if ~isempty(ii)
   cver = attval{ii};
end

%---------------------------------------------------------------------------
