function dpth = dat_path(vers)

global CSL_Version

if nargin==0 | isempty(vers)
   if isempty(CSL_Version)
      vers = 2;
   else
      vers = CSL_Version;
   end
end

switch vers
  case {1,2}
    dpth = path_pc_or_nix('eez_data/');
  case 3
    dpth = path_pc_or_nix('eez_data/boa_csl3/');
  otherwise
    dpth = '';
end
