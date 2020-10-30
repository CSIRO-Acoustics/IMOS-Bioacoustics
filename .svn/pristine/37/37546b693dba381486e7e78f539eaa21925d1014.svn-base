function nmea = parsenmea(text)
% parsenmea Parses a NMEA text string and returns a data structure with the
% values of interest returned.
%
% This only supports some sentence types.
%
% Gordon.Keith@csiro.au 20130605

if length(text) < 8 || text(1) ~= '$' || text(7) ~= ','
    if ~strcmp(text(1:6), '$PSXN,')                          % don't warn for PSXN whatever it is.
        warning ('Invalid NMEA sentence: "%s"', text);
    end
end

txt = strtrim(text);
txt(txt == 0) = [];
if txt(end-2) == '*'
    txt(end-2:end) = '';        % drop checksum
end
fields = regexp(txt,',', 'split');

nmea.talker = text(2:3);
nmea.sentence = text(4:6);

try
switch nmea.sentence
    case 'GLL'
        nmea.latitude  = todegrees(fields{2}, fields{3});
        nmea.longitude = todegrees(fields{4}, fields{5});
    case 'GGA'
        nmea.time = totime(fields{2});
        nmea.latitude  = todegrees(fields{3}, fields{4});
        nmea.longitude = todegrees(fields{5}, fields{6});
        
    case 'RMC'
        nmea.time = totime(fields{2}) + todate(fields{10});
        nmea.latitude  = todegrees(fields{4}, fields{5});
        nmea.longitude = todegrees(fields{6}, fields{7});
        
    case 'ZDA'
        nmea.time = datenum(str2double(fields{5}),str2double(fields{4}),str2double(fields{3})) + ...
            totime(fields{2});
        
    case 'RPY'
        nmea.roll  = str2double(fields{2});
        nmea.pitch = str2double(fields{3});
    %   nmea.yaw   = str2double(fields{4});
    
    case 'MOT'
        switch nmea.talker
            case 'PA' % CSIRO
                nmea.roll  = str2double(fields{3});
                nmea.pitch = str2double(fields{5});
            case 'YX' % SeaLord
                nmea.roll  = str2double(fields{4});
                nmea.pitch = str2double(fields{2});
                nmea.yaw   = str2double(fields{6});
        end
        
    case 'GIM'
        nmea.gimblex = str2double(fields{2});
        nmea.gimbley = str2double(fields{4});
        
    case 'DPT'
        nmea.depth=str2double(fields{2});
        
    case 'DPH'
        dph = fields{2};
        ast = find(dph == '*',1);
        if ~isempty(ast)
            dph = dph(1:ast-1);
        end
        nmea.depth = str2double(dph);
        
    case 'CTD'
        nmea.CTDtemp = str2double(fields{2});
        nmea.CTDcond = str2double(fields{4});
        nmea.CTDpres = str2double(fields{6});
        nmea.CTDsalt = str2double(fields{8});    
        if length(fields) > 10
            nmea.CTDoxy = str2double(fields{10});   
        end
        
end
catch e
    warning('Problem parsing NMEA sentence: \n%s\n"%s"', e.message, text);
end

function degrees = todegrees(coord, hemisphere)

coord = str2double(coord);
degrees = floor(coord / 100);
minutes = coord - 100 * degrees;
degrees = degrees + minutes / 60;
if hemisphere(1) == 'S' || hemisphere(1) == 'W'
    degrees = - degrees;
end

function time = totime(field)

hour = str2double(field(1:2));
min  = str2double(field(3:4));
sec  = str2double(field(5:end));
time = (hour + (min + sec / 60) / 60) / 24;

function date = todate(field)
if field(5) >= '7' && field(5) <= '9'
    cent = 1900;
else
    cent = 2000;
end
date = datenum(cent + str2double(field(5:6)),str2double(field(3:4)),str2double(field(1:2)));

