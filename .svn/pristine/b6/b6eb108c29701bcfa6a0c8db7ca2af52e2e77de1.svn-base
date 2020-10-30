% DECODE_GREYLIST  Prepare a user-friendly mat-file from the greylist text
% file downloaded from the Argo GDAC.
%
% Jeff Dunn
%
% NOTE:  Compare Argo serial_date with "rdate" by subtracting Matlab date
%        for 1/1/1900  
%       >> t1900 = datenum([1900 1 1 0 0 0]);
%       >> if (prof.serial_date-t1900 > rdate(ii,1)) && ...
%	       (prof.serial_date-t1900 < rdate(ii,2))
%
% USAGE: decode_greylist 

disp('The latest greylist should be downloaded from ')
disp('http://www.usgodae.org/ftp/outgoing/argo/ar_greylist.txt');
disp('and copied to /home/datalib/observations/argo/');
disp(' ')
disp('If this script does not complete successfully you may need to find')
disp('and remove problem lines from the downloaded file.')

fid=fopen('/home/datalib/observations/argo/ar_greylist.txt','r');

tline = fgetl(fid);   % Read off header line

ii = 0;
tline = fgetl(fid);
while ischar(tline)   
   C = textscan(tline,'%d%s%d%d%d','delimiter',',');
   ii = ii+1;
   wmoid(ii) = double(C{1}(1));
   jj = strmatch(deblank(upper(C{2}(1))),['PRES';'TEMP';'PSAL';'COND'],'exact');
   if isempty(jj)
      disp(['Unknown parameter: ' C{2}(1)])
   else
      parnm(ii) = jj;
   end
   for ll = 1:2
      if ~isempty(C{2+ll}) && C{2+ll}(1)~=0
	 yr = floor(C{2+ll}(1)/10000);
	 cut = rem(C{2+ll}(1),10000);
	 mo = floor(cut/100);
	 sec = rem(cut,100);
	 rdate(ii,ll) = greg2time(double([yr mo sec 0 0 0]));
      else
	 rdate(ii,ll) = inf;
      end
   end
   if ~isempty(C{5})
      flag(ii) = double(C{5}(1));
   else
      flag(ii) = 3;
   end
   tline = fgetl(fid);
end

fclose(fid);

note = {'wmoid: WMO id of float', ...
	'flag: standard Argo QC flag. Generally ignore all data on this list',...
	'rdate: Start and end of dubious data, decimal days since 1900',...
	'parnm: parameter code: 1=pressure 2=temperate 3=salinity 4=conductivity',...
	' ',['This file updated ' date]};

save /home/datalib/observations/argo/greylist wmoid flag rdate parnm note
