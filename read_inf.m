function meta = read_inf(meta, inf_file)
% Read an .inf file for the voyage and return voyage metadata.
%
% Inputs:
%   meta        structure containing metadata, may be empty
%   inf_file    full path to voyage information file (.inf for whole voyage)
%
% Outputs:
%   meta        structure containing metadata
%
% This function uses ports.txt to locate the nearest port in the BODC
% database.
%

start = 0;
fin = 0;

try
    % read inf file
    fid    = fopen(inf_file, 'rt');
    
    line = fgetl(fid);
    while (ischar(line))
        if length(line) > 4
            switch line(1:4)
                case 'Star'
                    start = 1;
                case 'End '
                    fin = 1;
                case 'Time'
                    if fin
                        meta.cruise_end_date = sscanf(line, 'Time: %s%c%s JD');
                    elseif start
                        meta.cruise_start_date = sscanf(line, 'Time: %s%c%s JD');
                    end
                case 'Lon:'
                    if fin
                        endpos =  sscanf(line, 'Lon: %f Lat: %f');
                    elseif start
                        startpos = sscanf(line, 'Lon: %f Lat: %f');
                    end
                case 'Mini'
                    lon = sscanf(line,'Minimum Longitude:%f Maximum Longitude: %f');
                    lat = sscanf(line,'Minimum Latitude:%f Maximum Latitude: %f');
                    if ~isempty(lon)
                        meta.cruise_westlimit = lon(1);
                        meta.cruise_eastlimit = lon(2);
                        meta.cruise_units = 'signed decimal degrees';
                        meta.cruise_projection = 'Geographic';
                    end
                    if ~isempty(lat)
                        meta.cruise_southlimit = lat(1);
                        meta.cruise_northlimit = lat(2);
                    end
            end
        end
        line = fgetl(fid);
    end
    
    fclose(fid);
catch e
    warning('read_inf:exception', 'Problem reading inf file %s\n%s', ...
        inf_file, e.message)
end

%
% Attempt to find port near start and end positions.
%
% ports.txt contains BODC codes in the format downloaded from 
% http://seadatanet.maris2.nl/v_bodc_vocab/browse_export.asp?l=C38&all=yes&description=1&x=39&y=12
% These have been restricted to valid southern hemisphere ports by the
% command (710 lines down from 4901):
% grep "latitude>-" vocab_C38.csv > ports.txt
%
% Further editing of the ports.txt file may be desirable (e.g. BSH1325
% Sandgate is closer to where the Surveyor ties up than BSH1324 Brisbane,
% it may be desirable to drop Sandgate).

try
    portfile = fullfile(fileparts(which('read_inf')), 'ports.txt');
    if exist(portfile, 'file') == 2
        pfid = fopen(portfile, 'rt');
        ports = textscan(pfid, '%q%q%*q"<country>%*s/country><latitude>%f/latitude><longitude>%f/longitude>%*[^\n]', ...
            'delimiter',';<','CommentStyle','#');
        fclose(pfid);
        
        [BODC,name]=closest_port(ports,startpos,'cruise start');
        if ~isempty(BODC)
            meta.cruise_start_BODC_code = BODC;
            meta.cruise_start_port = name;
        end
        [BODC,name]=closest_port(ports,endpos,'cruise end');
        if ~isempty(BODC)
            meta.cruise_end_BODC_code = BODC;
            meta.cruise_end_port = name;
        end
    end
catch e
    warning('read_inf:ports', 'Problem reading ports file %s\n%s', ...
        portfile, e.message)
end

function [code, name] = closest_port(ports, pos, prompt)

BODC=1; NAME=2; LAT=3; LON=4;

close = .05;    % 25 km
range = .2;     % 50 km

code='';
name='';

d2 = (ports{LAT} - pos(2)) .^2 + (ports{LON} - pos(1)) .^2;

[dist2, order] = sort(d2);

if dist2(1) < close && dist2(2) > range     % unique good candidate
    code = ports{BODC}{order(1)};
    name = ports{NAME}{order(1)};
else
    toofar = find(dist2 > range, 1);
    while toofar == 1
        range = range * 2;
        toofar = find(dist2 > range, 1);
    end
    near_ports = ports{NAME}(order(1:toofar));
    near_ports{end} = 'more ...';
    ok = 1;
    while ok
        [selection, ok] = listdlg('ListString', near_ports, ...
            'SelectionMode', 'single', ...
            'InitialValue', 1, ...
            'Name', 'Port Selection', ...
            'PromptString', ['Please select ' prompt ' port'], ...
            'ListSize', [400 300]);
        if ok && selection < length(near_ports)    % more
            port = order(selection);
            code = ports{BODC}{port};
            name = ports{NAME}{port};
            ok = 0;
        end
        if ok
            range = range * 2;
            toofar = find(dist2 > range, 1);
            near_ports = ports{NAME}(order(1:toofar));
            near_ports{end} = 'more ...';
        end
    end
        
end
