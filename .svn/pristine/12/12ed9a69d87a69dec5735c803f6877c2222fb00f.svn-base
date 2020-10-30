% Reading EK60 raw data file configuration header
% Simrad, Lars Nonboe Andersen, 21/12-01
% modified Gordon.Keith@csiro.au 20130528

function configheader = readconfigheader(fid)

configheader.surveyname         = fread(fid,128,'*char')';
configheader.transectname       = fread(fid,128,'*char')';
configheader.soundername        = fread(fid,128,'*char')';

configheader.spare              = fread(fid,128,'*char')';

configheader.transducercount    = fread(fid,1,'int32');

