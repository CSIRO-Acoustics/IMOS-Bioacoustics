% SST_FILNAM  Append date-appropriate postfix to an SST filename.
%
%   Some CSIRO SST files have yearly files up to 2004, then 6 monthly
%   files after that. This gets trickier if we want the file before or
%   after the one containing the given time. 
%
% INPUT
%  tim     Required time, in days-since-1/1/1900
%  fnm     Base filename
%  nxt     +1 or -1 if want the next or previous file [default 0]
%
% OUTPUT
%  fnmo    The extended filename
%
% Jeff Dunn  CMAR  15/8/2007
%
% CALLED BY:  get_sst.m   get_sst_xy.m  get_sst_xyt.m
%
% USAGE: fnmo = sst_filnam(tim,fnm,nxt);

function fnmo = sst_filnam(tim,fnm,nxt)

greg = time2greg(tim);
yr = greg(1); 
mon = greg(2);

if nargin==3 && ~isempty(nxt)
   if nxt==1
      if yr<=2004 || mon>=7
	 yr = yr+1;
	 mon = 1;
      else
	 mon = 7;
      end
   elseif nxt==-1
      if yr<=2004 || mon<=6
	 yr = yr-1;
	 mon = 12;
      else
	 mon = 6;
      end
   end
end

if yr<=2004
   fnmo = [fnm num2str(yr)];
elseif mon<=6
   fnmo = [fnm num2str(yr) '01-06'];
else
   fnmo = [fnm num2str(yr) '07-12'];
end	       

%--------------------------------------------------------------------------
