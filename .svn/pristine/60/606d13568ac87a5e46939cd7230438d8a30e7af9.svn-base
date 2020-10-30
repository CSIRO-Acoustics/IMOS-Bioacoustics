% DSET_DUP_PREF  return default preferred datasets for duplicate removal, in
%    the required cell-array for GET_ALL_CSL
%
% INPUT:  
%   dsets  vector of dataset codes   (see DSET_NAME for meanings of codes)
%   var    the property code  
%   dset1  dataset code - the single dataset whose casts are being thinned
%          (if omit then dupref has a cell for each dataset in dsets.)
% OUTPUT:  
%  if no dset1, cell array, one cell per dataset code, containing vector of
%          codes of datasets whose duplicates are preferred to this dataset.
%  if dset1, dupref is just a vector.
%
% Jeff Dunn  15/10/04 
%
% USAGE: dupref = dset_dup_pref(dsets,var,dset1);

function dupref = dset_dup_pref(dsets,var,dset1)

% TEST METHOD: 
%   To ensure that 2 datasets are not both excluded because of the other:
%   
%   allds = [7:14 17:29 31 34 35 70:72 81 121];
%   zz = zeros(max(allds));   
%   for ii = allds
%      pp = dset_dup_pref(allds,1,ii);
%      if ~isempty(pp)
%         zz(ii,pp) = 1;
%      end
%   end 
%
%   [i1,i2] = find(zz==1 & zz'==1);
%
% Many datasets have no overlap, hence no need to exclude or be excluded by, 
% other datasets, and there is no efficient test for where datasets wrongly
% allow each other.

if nargin<3
   dset1 = [];
end

% Note: an empty cell for a dataset means it is a primary dataset - its data
%  is preferred to all others (or it could mean there are no other datasets 
%  which could have copies of the data in this one.) 

% NOTE: still need to work out prefs for WOD01 minor subsets (25-28)

% The listed datasets in defs{N} are those used in preference to dset N
% So, where the list is created using "setxor(allds,[N1 N2 ...]) then
% datasets N1 N2 are the only ones that are LESS favoured.

allds = [7:14 17:29 31 34 35 36 41:49 51:59 70:72 81 121];

defs{199} = [];

%if var==1
if var<3
   % T primary datasets 11 12 13 31 44 46
   defs{7}  = 19;
   defs{8}  = 11;
   defs{9}  = 19;
   defs{10} = setxor(allds,10);
   defs{14} = setxor(allds,[14 19 20]);
   defs{17} = [7 8 9 11 12 13 21:27 31];
   defs{18} = [7 8 9 11 12 13 17 19 20 21:28];
   defs{19} = [13 21:28 41 42];
   defs{20} = [7 8 9 11 12 13 19 21:28];
   defs{21} = [7 9 11 13 41];
   defs{22} = [7 9 11 13 21 41 42 72];
   defs{23} = 11;
   defs{24} = [7 9 11 13];
   defs{28} = [8 17];
   defs{35} = [7 9 13 19 21:24];
   defs{38} = [7 8 9 11 12 13 17 21:27 31];
   defs{41} = [7 9 11 13];
   defs{42} = [7 9 11 13 41 72];
   defs{43} = 11;
   defs{51} = [7 9 11 13];
   defs{52} = [7 9 11 13 51 72];
   defs{53} = 11;
   defs{71} = [7 9 13 19 21:24 41 42];
   defs{72} = [7 9 13 19 21 41];
elseif var==2
   % S primary datasets 11 12 13 31 44 46
   defs{7}  = [19 41];
   defs{9}  = [19 41];
   defs{10} = [7 9 13 19 21:27];
   defs{19} = [13 21:28 41:46];
   defs{21} = [7 9 11 13 41 42];
   defs{22} = [7 9 11 13 21 41 42];
   defs{23} = 11;
   defs{24} = [7 9 11 13];
   defs{41} = [11 13];
   defs{42} = [7 9 11 13];
   defs{43} = 11;
   defs{51} = [11 13];
   defs{52} = [7 9 11 13 51];
   defs{53} = 11;
   defs{35} = [7 9 13 19 21:24];
   defs{71} = [7 9 13 19 21:24];
   defs{72} = [7 9 13 19 21];
elseif var==3
   defs{10} = [7 9 13 19 21:27 51:54];
   defs{19} = [7 9 13 21:28 51:54];
   defs{21} = [7 9 13 22 41 42 51 52];
   defs{22} = [7 9 13 41 42 51 52];
   defs{24} = [7 9 13 41 42 51 52];
   defs{41} = [7 9 13 51];
   defs{42} = [7 9 13 41 51];
   defs{51} = [7 9 13];
   defs{52} = [7 9 13 51];
   defs{72} = [7 9 13 19 21 41 51 52];
else
   defs{10} = [7 9 13 19 21:27 52];
   defs{19} = [7 9 13 21:28 52];
   defs{21} = [7 9 13 41 52];
   defs{24} = [7 9 13 41 52];
   defs{41} = [7 9 13];
   if var==6
      % No longer automatically include 71 when extracting dset 7, so now
      % need to explicitly remove dset 71 dups from dset 22.
      defs{22} = [7 9 13 42 52];
      defs{42} = [7 9 13];
      defs{52} = [7 9 13];
      defs{71} = [19 22];
   else
      defs{22} = [7 9 13 42 52];      
      defs{42} = [7 9 13];      
      defs{52} = [7 9 13];      
   end
   defs{72} = [7 9 13 19];
end
   
nds = length(dsets);

if isempty(dset1)
   for ii = 1:nds
      kk = ismember(dsets,defs{dsets(ii)});
      dupref{ii} = dsets(kk);
   end
else
   kk = ismember(dsets,defs{dset1});
   dupref = dsets(kk);
end

%-------------------------------------------------------------------------
