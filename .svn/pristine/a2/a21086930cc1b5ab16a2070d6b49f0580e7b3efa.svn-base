% writing EK60 raw data file datagram header
% based on readdgheader - Simrad, Lars Nonboe Andersen, 8/5-03

function writedgheader(fid2,dgheader);

%dgheader.datagramtype = char(fread(fid,4,'char')');
%fwrite(fid2,cellstr(dgheader.datagramtype),'char'); % convert to cell str and write
fwrite(fid2,(dgheader.datagramtype),'char'); %  

%lowdatetime = fread(fid,1,'uint32');
fwrite(fid2,dgheader.lowdatetime,'uint32');

%highdatetime = fread(fid,1,'uint32');
fwrite(fid2,dgheader.highdatetime,'uint32');

%dgheader.datetime = highdatetime*2^32 + lowdatetime;
