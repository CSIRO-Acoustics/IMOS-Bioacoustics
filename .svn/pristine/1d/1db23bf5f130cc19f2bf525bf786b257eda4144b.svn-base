% READ_CSL_ARGO: Instead of using GET_ALL_CSL3 to extract some Argo data,
%    this is simpler if just want to read the whole CSL Argo dataset, but still
%    apply BOA screening.
%
% WARNING: This is a script, and will create and modify variables in your
% workspace. Safest to just copy this code and include or modify to suit your
% needs. The variables created or clobbered are:  
%   fnmi 
%   II
%   infodir
%   jdep
%  From Argo file: botdep cru lat lon mld prno stnno s s_castflag t t_castflag
%                  time mirror_date
%
% NOTE:  Where whole profiles are flagged as bad they will not be removed but
%        just filled with NaNs
%
% Jeff Dunn 21/3/2012

fnmi = '/home/eez_data/boa_csl3/argo_csl3';
load(fnmi,'botdep','cru','lat','lon','t','s', ...
     'mld','s_castflag','stnno','t_castflag',...
     'time','mirror_date');

% Apply all local screening
infodir = '/home/eez_data/boa_qc_data/';
load([infodir '11/t_scr'],'scrstn');

t(t_castflag>0,:) = nan; 
for jdep = 1:size(t,2)
   II = find(ismember(stnno,scrstn{jdep}));
   if ~isempty(II)
      t(II,jdep) = nan;
   end
end

s(s_castflag>0,:) = nan; 
load([infodir '11/s_scr'],'scrstn');
for jdep = 1:size(s,2)
   II = find(ismember(stnno,scrstn{jdep}));
   if ~isempty(II)
      s(II,jdep) = nan;
   end
end

