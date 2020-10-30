% Reading EK60 raw data file configuration transducer
% Simrad, Lars Nonboe Andersen, 8/5-03
% modified Gordon.Keith@csiro.au 20130528

function configtransducer = readconfigtransducer(fid)
configtransducer.channelid                      = fread(fid,128,'*char')';
configtransducer.beamtype                       = fread(fid,1,'int32');
configtransducer.frequency                      = fread(fid,1,'float32');
configtransducer.gain                           = fread(fid,1,'float32');
configtransducer.equivalentbeamangle            = fread(fid,1,'float32');
configtransducer.beamwidthalongship             = fread(fid,1,'float32');
configtransducer.beamwidthathwartship           = fread(fid,1,'float32');
configtransducer.anglesensitivityalongship      = fread(fid,1,'float32');
configtransducer.anglesensitivityathwartship    = fread(fid,1,'float32');
configtransducer.anglesoffsetalongship          = fread(fid,1,'float32');
configtransducer.angleoffsetathwartship         = fread(fid,1,'float32');
configtransducer.posx                           = fread(fid,1,'float32');
configtransducer.posy                           = fread(fid,1,'float32');
configtransducer.posz                           = fread(fid,1,'float32');
configtransducer.dirx                           = fread(fid,1,'float32');
configtransducer.diry                           = fread(fid,1,'float32');
configtransducer.dirz                           = fread(fid,1,'float32');
configtransducer.pulselengthtable               = fread(fid,5,'float32');
configtransducer.spare2                         = fread(fid,8,'*char')';
configtransducer.gaintable                      = fread(fid,5,'float32');
configtransducer.spare3                         = fread(fid,8,'*char')';
configtransducer.sacorrectiontable              = fread(fid,5,'float32');
configtransducer.spare4                         = fread(fid,52,'*char')';
