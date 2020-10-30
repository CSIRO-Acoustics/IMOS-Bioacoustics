% DSET_NAME  Return dataset name(s) corresponding to dataset codes (or if
%  no argument given, list all codes and names).
%
% INPUT: dsets:  dataset numeric codes
%        shrt:   0=full names,  1=short names  [def 0]
%
% USAGE: dnms = dset_name(dsets,shrt);

function dnms = dset_name(dsets,shrt)

fnm = path_pc_or_nix('eez_data/software/matlab/dset_names');
load(fnm)

if nargin<2 | isempty(shrt)
   shrt = 0;
end

if nargin<1 | isempty(dsets)
   for ii = 1:length(dset_names)
      if ~isempty(dset_names{ii})
	 if shrt
	    disp([num2str(ii)  '   ' dset_short_names{ii}]);
	 else
	    disp([num2str(ii)  '   ' dset_names{ii}]);
	 end
      end
   end
   dnms = [];
else
   if shrt
      dnms = dset_short_names(dsets);
   else
      dnms = dset_names(dsets);
   end
end
