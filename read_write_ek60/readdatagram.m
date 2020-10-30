function datagram = readdatagram( fid, varargin)
%READDATAGRAM reads an EK60 datagram for the file ID fid.
%
% Inputs:
%   fid     ID of an open EK60 file
%   'noraw'     Don't read RAW0 datagrams
%   'rawheader' Only read the header of RAW0 datagrams
%   'fullraw'   Fully read RAW0 datagrams (default)
%   'nonmea'    Don't parse NME0 datagrams
%   'nmea'      Parse NME0 datagrams into a structure (default)
%
% Output
%   datagram    structure containing datagram contents, empty at end of
%               file

% read controls from varargin

NONE    = 0;
HEADER  = 1;
FULL    = 2;

HEADERLENGTH = 12;

detail = FULL;
parse = 1;

for i = 1:nargin - 1
    switch lower(varargin{i})
        case 'noraw'
            detail = NONE;
            
        case 'rawheader'
            detail = HEADER;
            
        case 'fullraw'
            detail = FULL;
            
        case 'nonmea'
            parse = 0;
            
        case 'nmea'
            parse = 1;
    end
end

% read length of datagram

len = fread(fid,1,'uint32');

if isempty(len) % EOF
    datagram = [];
    return
end
    
% read datagram header to determine type
datagram = readdgheader(fid);

% read datagram by type
switch datagram.datagramtype
    case 'CON0'
        datagram.CON0 = readconfigheader(fid);
        
        for i = 1:datagram.CON0.transducercount
            datagram.CON0.transducer(i) = readconfigtransducer(fid);
        end
        
    case 'CON1'
	datagram.text = fread(fid,len - HEADERLENGTH, '*char')';

    case 'NME0'
        datagram.text = fread(fid,len - HEADERLENGTH, '*char')';
        if parse
            datagram.NME0 = parsenmea(datagram.text);
        end
    case 'TAG0'
        datagram.text = fread(fid,len - HEADERLENGTH, '*char')';
        
    case 'RAW0'
        if nargin < 2
            detail = FULL;
        end
        if detail == NONE
            fseek(fid, len - HEADERLENGTH, 0);
        else                
            datagram.RAW0 = readsampledata(fid, detail ~= HEADER);
        end
        
    otherwise
        error(['Unknown datagram ' datagram.datagramtype]);
end

if len ~= fread(fid,1,'uint32')
    error('Corrupt file, record length mismatch');
end


