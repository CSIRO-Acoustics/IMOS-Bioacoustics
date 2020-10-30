function nmea = readnmea(fid,len, parse)

text = fread(fid,len, '*char')';
nmea.text = text;

if length(text) < 8 || text(1) ~= '$' || text(7) ~= ','
    error ([ 'Invalid NMEA sentence ' text]);
end

if nargin < 3 || parse
    
    fields = regexp(text,',', 'split');
    
    nmea.talker = text(2:3);
    nmea.sentence = text(4:6);
    
    switch nmea.sentence
        case 'GLL'
            nmea.latitude  = todegrees(fields(2), fields(3));
            nmea.longitude = todegrees(fields(4), fields(5));
        case 'GGA'
            nmea.time = totime(fields(2));
            nmea.latitude  = todegrees(fields(3), fields(4));
            nmea.longitude = todegrees(fields(5), fields(6));
            
        case 'RMC'
            nmea.time = totime(fields(2)) + todate(fields());
            nmea.latitude  = todegrees(fields(4), fields(5));
            nmea.longitude = todegrees(fields(6), fields(7));
            
        case 'ZDA'
            nmea.time = datenum(fields(3),fields(4),fields(5)) + ...
                totime(fields(2));
            
        case 'RPY'
            nmea.roll  = str2double(fields(2));
            nmea.pitch = str2double(fields(3));
        %   nmea.yaw   = str2double(fields(4));
    end
end


function degrees = todegrees(coord, hemisphere)

coord = str2double(coord);
degrees = floor(coord / 100);
minutes = coord - 100 * degrees;
degrees = degrees + minutes / 60;
if hemisphere == 'S' || hemisphere == 'W'
    degrees = - degrees;
end

function time = totime(field)

hour = str2double(field(1:2));
min  = str2double(field(3:4));
sec  = str2double(field(5:end));
time = (hour + (min + sec / 60) / 60) / 24;

function date = todate(field)
date = datenum(2000 + field(1:2),field(3:4),field(5:6);