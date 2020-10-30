% writing EK60 raw data file configuration transducer
% ADR - modfied from readconfigtransducer

function writeconfigtransducerraw(fid2,configtransducer);

%configtransducer.channelid = char(fread(fid,128,'char')');
fwrite(fid2,configtransducer.channelid,'char'); % raw version
%fwrite(fid2,configtransducer.channelid,'int32'); % raw version
%configtransducer.beamtype = fread(fid,1,'int32');
fwrite(fid2,configtransducer.beamtype,'int32');
%configtransducer.frequency = fread(fid,1,'float32');
fwrite(fid2,configtransducer.frequency,'float32');
%configtransducer.gain = fread(fid,1,'float32');
fwrite(fid2,configtransducer.gain,'float32');
%configtransducer.equivalentbeamangle = fread(fid,1,'float32');
fwrite(fid2,configtransducer.equivalentbeamangle,'float32');
%configtransducer.beamwidthalongship = fread(fid,1,'float32');
fwrite(fid2,configtransducer.beamwidthalongship,'float32');
%configtransducer.beamwidthathwartship = fread(fid,1,'float32');
fwrite(fid2,configtransducer.beamwidthathwartship,'float32');
%configtransducer.anglesensitivityalongship = fread(fid,1,'float32');
fwrite(fid2,configtransducer.beamwidthathwartship,'float32');
%configtransducer.anglesensitivityathwartship = fread(fid,1,'float32');
fwrite(fid2,configtransducer.anglesensitivityathwartship,'float32');
%configtransducer.anglesoffsetalongship = fread(fid,1,'float32');
fwrite(fid2,configtransducer.anglesoffsetalongship,'float32');
%configtransducer.angleoffsetathwartship = fread(fid,1,'float32');
fwrite(fid2,configtransducer.angleoffsetathwartship,'float32');
%configtransducer.posx = fread(fid,1,'float32');
fwrite(fid2,configtransducer.posx,'float32');
%configtransducer.posy = fread(fid,1,'float32');
fwrite(fid2,configtransducer.posy,'float32');
%configtransducer.posz = fread(fid,1,'float32');
fwrite(fid2,configtransducer.posz,'float32');
%configtransducer.dirx = fread(fid,1,'float32');
fwrite(fid2,configtransducer.dirx,'float32');
%configtransducer.diry = fread(fid,1,'float32');
fwrite(fid2,configtransducer.diry,'float32');
%configtransducer.dirz = fread(fid,1,'float32');
fwrite(fid2,configtransducer.dirz,'float32');
%configtransducer.pulselengthtable = fread(fid,5,'float32');
fwrite(fid2,configtransducer.pulselengthtable,'float32');
%configtransducer.spare2 = char(fread(fid,8,'char')');
fwrite(fid2,configtransducer.spare2,'char'); % raw version
%configtransducer.gaintable = fread(fid,5,'float32');
fwrite(fid2,configtransducer.gaintable,'float32')
%configtransducer.spare3 = char(fread(fid,8,'char')');
fwrite(fid2,configtransducer.spare3,'char');% raw version 
%configtransducer.sacorrectiontable = fread(fid,5,'float32');
fwrite(fid2,configtransducer.sacorrectiontable,'float32')
%configtransducer.spare4 = char(fread(fid,52,'char')');
fwrite(fid2,configtransducer.spare4,'char') % raw version