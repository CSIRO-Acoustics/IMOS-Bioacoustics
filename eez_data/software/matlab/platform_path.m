% PLATFORM_PATH
% Construct paths appropriately for the type of machine being used.
%
%  **** SUPERCEDED BY        nix_px_path.m            ****
%
% INPUT: mach  - home of the disc, so can construct PC-style path
%        upth  - unix style path, without the '/home/' 
%
% OUTPUT  pth  - complete path
%        slsh  - slash symbol, if want to construct extra paths
%        plat  - 1=PC, 0=Unix, -1=MAC
%
% EG:  pth = platform_path('reg','reg2/SST_mw/netcdf/');
%
% USAGE: [pth,slsh,plat] = platform_path(mach,upth);

function [pth,slsh,plat] = platform_path(mach,upth)

% MODS: 21/3/2012 Mucked around because new Unix paths don't often start with
%          "/home" so have to allow whole path (as seen by a PC) to be specified.

cname = computer;
if strncmp(cname,'PC',2)
   plat = 1;
   if isempty(strfind(mach,'-'))
      pth = ['\\' mach '-hf\'];
   else
      pth = ['\\' mach '\'];
   end      
   ii = findstr(upth,'/');
   slsh = '\';
   upth(ii) = slsh;
   pth = [pth upth];
elseif strncmp(cname,'MAC',3)
   plat = -1;
   disp([7 'Sorry - do not how to find datafiles from a Mac'])
   pth = '';
   slsh = '?';
else
   % Assuming not a VAX, must be Unix
   plat = 0;
   if upth(1)=='/'       
      pth = upth;
   else
      pth = ['/home/' upth];
   end
   slsh = '/';
end


return

%-----------------------------------------------------------------------------
