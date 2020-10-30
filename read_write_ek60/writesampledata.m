% writing raw data file sample data
% mod ADR to write raw data
function writesampledataraw(fid2,sampledata)
%sampledata.channel = fread(fid,1,'int16');
fwrite(fid2,sampledata.channel,'int16');
%sampledata.mode_low = fread(fid,1,'int8');
fwrite(fid2,sampledata.mode_low,'int8');
%sampledata.mode_high = fread(fid,1,'int8');
fwrite(fid2,sampledata.mode_high,'int8');
%sampledata.mode = 256*sampledata.mode_high + sampledata.mode_low;
%sampledata.transducerdepth = fread(fid,1,'float32');
fwrite(fid2,sampledata.transducerdepth,'float32');
%sampledata.frequency = fread(fid,1,'float32');
fwrite(fid2,sampledata.frequency,'float32');
%sampledata.transmitpower = fread(fid,1,'float32');
fwrite(fid2,sampledata.transmitpower,'float32');
%sampledata.pulselength = fread(fid,1,'float32');
fwrite(fid2,sampledata.pulselength,'float32');
%sampledata.bandwidth = fread(fid,1,'float32');
fwrite(fid2,sampledata.bandwidth,'float32');
%sampledata.sampleinterval = fread(fid,1,'float32');
fwrite(fid2,sampledata.sampleinterval,'float32');
%sampledata.soundvelocity = fread(fid,1,'float32');
fwrite(fid2,sampledata.soundvelocity,'float32');
%sampledata.absorptioncoefficient = fread(fid,1,'float32');
fwrite(fid2,sampledata.absorptioncoefficient,'float32');
%sampledata.heave = fread(fid,1,'float32');
fwrite(fid2,sampledata.heave,'float32');
%sampledata.roll = fread(fid,1,'float32');
fwrite(fid2,sampledata.roll,'float32');
%sampledata.pitch = fread(fid,1,'float32');
fwrite(fid2,sampledata.pitch,'float32');
%sampledata.temperature = fread(fid,1,'float32');
fwrite(fid2,sampledata.temperature,'float32');
%sampledata.trawlupperdepthvalid = fread(fid,1,'int16');
fwrite(fid2,sampledata.trawlupperdepthvalid,'int16');
%sampledata.trawlopeningvalid = fread(fid,1,'int16');
fwrite(fid2,sampledata.trawlopeningvalid,'int16');
%sampledata.trawlupperdepth = fread(fid,1,'float32');
fwrite(fid2,sampledata.trawlupperdepth,'float32');
%sampledata.trawlopening = fread(fid,1,'float32');
fwrite(fid2,sampledata.trawlopening,'float32');
%sampledata.offset = fread(fid,1,'int32');
fwrite(fid2,sampledata.offset,'int32');
%sampledata.count = fread(fid,1,'int32');
fwrite(fid2,sampledata.count,'int32');
%sampledata.raw_power = fread(fid,sampledata.count,'int16');
fwrite(fid2,sampledata.raw_power,'int16');

%sampledata.power = sampledata.raw_power*10*log10(2)/256;  % converted data
if (sampledata.mode>1)
    %sampledata.raw_angle = fread(fid,[2 sampledata.count],'int8');  % this is probably my problem
    fwrite(fid,sampledata.raw_angle,'int8');
    %angle=sampledata.raw_angle;
    %sampledata.angle = angle(1,:) + angle(2,:)*256;
    %sampledata.alongship = angle(2,:)';
    %sampledata.athwartship = angle(1,:)';
end

