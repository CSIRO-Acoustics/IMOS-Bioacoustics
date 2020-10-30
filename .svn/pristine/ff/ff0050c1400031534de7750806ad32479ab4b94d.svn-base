% NIX_PC_PATH    
%
%   ***  Deprecated - use PATH_PC_OR_NIX instead
%
% Construct paths appropriately for the type of machine being used.
%
% INPUT: mach  - suffix to construct PC-style path (leave off leading and
%                trailing '\')
%        upth  - unix style path, without the '/home/'. Leave off leading '/'
%                unless path does not start with '/home/'
%
% OUTPUT  pth  - complete path
%        slsh  - slash symbol, if want to construct extra paths
%        plat  - 1=PC, 0=Unix, -1=MAC
%
% EG:  pth = nix_pc_path('fstas-hba\CMAR\Project1','reg2/SST_mw/netcdf/');
%
% DEVOLVED from and replaces platform_path 26/3/2012
%
% NOTE: If samba is not running on a machine (so can't see from PC) then you
%    can create a symbolic link in a samba visible directory (on a
%    samba-running  machine)  to the NFS directory you want then samba will 
%    follow the link.  (according to Gordon.)
%
% USAGE: [pth,slsh,plat] = nix_pc_path(mach,upth);

function [pth,slsh,plat] = nix_pc_path(mach,upth)

% Paths as at 26 Mar 2012
%  
% fstas2-hba:/CMAR/Share1/UOT-data
% fstas2-hba:/CMAR/Share/unix-apps/matlab/R2008b
% fstas2-hba:/CMAR-HOME3/dunn
% fstas2-hba:/datalib
% fstas2-hba:/CMAR/Project1/argo
% fstas2-hba:/CMAR/Project1/UOT
% fstas2-hba:/CMAR/Project1/netcdf-data
% fstas2-hba:/CMAR/Project1/eez_data
% oceania-hf:/work/oez5
%
% As at 25 May 2012, for PCs these paths may have \CSIRO in front of \CMAR ??

cname = computer;
if strncmp(cname,'PC',2)
   plat = 1;
   pth = ['\\' mach '\'];
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
