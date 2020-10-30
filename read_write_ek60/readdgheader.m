% Reading EK60 raw data file datagram header
% Simrad, Lars Nonboe Andersen, 8/5-03
% mod by ADR 2/17/03
% modified Gordon.Keith@csiro.au 20130528

function dgheader = readdgheader(fid)

dgheader.datagramtype = fread(fid,4,'*char')';

dgheader.lowdatetime  = fread(fid,1,'uint32');
dgheader.highdatetime = fread(fid,1,'uint32');

dgheader.datetime = (dgheader.highdatetime*2^32 + dgheader.lowdatetime) ...
    / 8.6400e+11 + datenum(1601,1,1);

