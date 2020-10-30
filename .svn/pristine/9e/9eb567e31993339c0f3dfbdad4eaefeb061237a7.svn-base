% SHOW_DUP_PREFS  Show the default preferred datasets for duplicate removal
%
% INPUT:  
%   dset   dataset code of the single dataset of interest
%   var    the property code (dups relations may differ for different variables)
%
% CALLS: dset_dup_pref
%
% Jeff Dunn  26/09/06 
%
% USAGE: show_dup_prefs(dset,var)

function show_dup_prefs(dset,var)

if var==1
   asets = [7 8 9 10 11 12 13 14 17 18 19 20 21:28 31 34 35 41:46 51:54];
elseif var==2
   asets = [7 9 10 11 12 13 19 21:27 31 35 41:46 51:54];
elseif var==3
   asets = [7 9 13 19 21:22 31 41 42 51 52];
elseif var==6
   asets = [7 9 13 19 22 31 71 72 42 52];
else
   asets = [7 9 13 19 22 31 72 42 52];
end

ii = find(asets==dset);
if isempty(ii)
   disp('Maybe this dataset does not contain this variable? Talk to Jeff')
else
   asets(ii) = [];
   dd = dset_dup_pref(asets,var,dset);
   if isempty(dd)
      disp('This is a "primary" dataset - it is preferred to all others')
   else
      disp('Profiles in this dataset will NOT be used if they occur in ')
      disp(['datasets ' int2str(dd)])
      if length(dd)<length(asets)
	 cc = setxor(asets,dd);
	 disp('Profiles in this dataset WILL be used in preference to those')
	 disp(['in datasets ' int2str(cc)])
      end	 
   end
end

%-------------------------------------------------------------------------
