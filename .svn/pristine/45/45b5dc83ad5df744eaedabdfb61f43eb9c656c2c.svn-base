% writing  EK60 raw data file configuration header
% ADR 2/14/04 based on readconfigheader 
% Simrad, Lars Nonboe Andersen, 21/12-01

function writeconfigheader(fid2,configheader);

%configheader.surveyname = char(fread(fid,128,'char')');
fwrite(fid2,(configheader.surveyname),'char'); % convert to cell str and write
%configheader.transectname = char(fread(fid,128,'char')');
fwrite(fid2,(configheader.transectname),'char'); % convert to cell str and write
%configheader.soundername = char(fread(fid,128,'char')');
fwrite(fid2,(configheader.soundername),'char'); % convert to cell str and write
%configheader.spare = char(fread(fid,128,'char')');
fwrite(fid2,(configheader.spare),'char'); % convert to cell str and write
%configheader.transducercount = fread(fid,1,'int32');
fwrite(fid2,configheader.transducercount,'int32'); % convert to cell str and write