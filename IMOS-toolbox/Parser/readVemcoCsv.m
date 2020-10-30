function [data, comment] = readVemcoCsv(dataLines, procHeader)
%readVemcoCsv Processes data section from a Vemco .csv file.
%
% This function is able to process data retrieved from a converted (.csv)
% data file generated by the Vemco Vue Logger program. This
% function is called from VemcoParse. Code modelled on readSBE37cnv.m
%
% Inputs:
%   dataLines - cell array of columns of raw data.
%   procHeader - Struct containing processed header.
%
% Outputs:
%   data       - Struct containing variable data.
%   comment    - Struct containing variable comment.
%
% Author: 		Simon Spagnol <s.spagnol@aims.gov.au>
% Contributor: 	Guillaume Galibert <guillaume.galibert@utas.edu.au>

%
% Copyright (C) 2017, Australian Ocean Data Network (AODN) and Integrated 
% Marine Observing System (IMOS).
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.
% If not, see <https://www.gnu.org/licenses/gpl-3.0.en.html>.
%
narginchk(2,2);

data = struct;
comment = struct;

columns = procHeader.columns;
% assume date and time would always be the first and second column
iDate=1;
iTime=2;
iDateTimeCol=[iDate iTime];
iProcCol=setdiff((1:length(columns)),iDateTimeCol);

% I don't know how to handle seperate date/time column in loop nicely
% so pull out datetime and set here, and process the other columns in the
% loop.
data.TIME = datenum(strcat(dataLines{iDate},'T',dataLines{iTime}),'yyyy-mm-ddTHH:MM:SS');
comment.TIME = 'TIME';

for kk = 1:length(iProcCol)
    iCol=iProcCol(kk);
    
    d = dataLines{iCol};
    
    [n, d, c] = convertData(genvarname(columns{iCol}), d);
    
    if isempty(n) || isempty(d), continue; end
    
    % if the same parameter appears multiple times,
    % don't overwrite it in the data struct - append
    % a number to the end of the variable name, as
    % per the IMOS convention
    count = 0;
    nn = n;
    while isfield(data, nn)
        
        count = count + 1;
        nn = [n '_' num2str(count)];
    end
    
    data.(nn) = d;
    comment.(nn) = c;
end

end

function [name, data, comment] = convertData(name, data)
%CONVERTDATA In order to future proof the .csv file, utilize the same ideal
% as for reading SBE37 data. This function is just a big switch statement which takes
% column header as input, and attempts to convert it to IMOS compliant name and
% unit of measurement. Returns empty string/vector if the parameter is not
% supported.

switch name
    
    %'Temperature (�C)'
    case {'Temperature0x280xFFFDC0x29', 'Temperature0x280xB0C0x29'};
        name = 'TEMP';
        comment = '';
        
    otherwise
        name = '';
        data = [];
        comment = '';
end
end