% TIME2DOY Convert time (as decimal days since the start of any year) to 
% day-of-year.
%
% J Dunn 28/3/96
%
% USAGE: doy = time2doy(time);

function doy = time2doy(time)

gtime = gregorian(time+2415020.5);

date0 = [NaN 1 1 0 0 0];
date0 = date0(ones([length(time) 1]),:);
date0(:,1) = gtime(:,1);
day0 = julian(date0);
doy = julian(gtime) - day0;

% % Just get doy the same shape as time.

[r c] = size(time);
if r==1; doy = doy'; end

% Paul Barker suggested that the code above is overkill, and algorithm below
% is sufficient and also base-year independant. However, Some folks *DO*
% want precise d-o-y, and that below is NOT precise! 15/2/08
%
% doy = time - floor(time/365.25)*365.25;


