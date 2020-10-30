function sample_data = timeDriftPP(sample_data, qcLevel, auto)
%TIMEDRIFTPP Prompts the user to apply time drift correction to the given data 
% sets. A pre-deployment time offset and end deployment time offset are
% required and if included in the DDB (or CSV file), will be shown in the
% dialog box. Otherwise, user is required to enter them.
%
% Offsets should be entered in seconds from UTC time (instrumentTime - UTCtime).
% An offset at the start will result in an offset to the start time. Global
% attributes of time coverage are also adjusted.
% Time drift calculations are applied to both raw and QC datasets.
% All IMOS datasets should be provided in UTC time. Raw data may not
% necessarily have been captured in UTC time, so a correction must be made
% before the data can be considered to be in an IMOS compatible format.
%
%
% Inputs:
%   sample_data - cell array of structs, the data sets to which time
%                 correction should be applied.
%   qcLevel     - string, 'raw' or 'qc'. Some pp not applied when 'raw'.
%   auto        - logical, run pre-processing in batch mode.
%
% Outputs:
%   sample_data - same as input, with time correction applied.
%

%
% Author:       Rebecca Cowley <rebecca.cowley@csiro.au>
% Contributor:  Guillaume Galibert <guillaume.galibert@utas.edu.au>
%

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

narginchk(2,3);

if ~iscell(sample_data), error('sample_data must be a cell array'); end
if isempty(sample_data), return;                                    end

% no modification of data is performed on the raw FV00 dataset except
% local time to UTC conversion
if strcmpi(qcLevel, 'raw'), return; end

% auto logical in input to enable running under batch processing
if nargin<3, auto=false; end

descs           = {};
nSample         = length(sample_data);
startOffsets    = zeros(nSample, 1);
endOffsets      = startOffsets;
sets            = ones(nSample, 1);

% create descriptions, and get timezones/offsets for each data set
for k = 1:nSample
    
    descs{k} = genSampleDataDesc(sample_data{k});
    
    %check to see if the offsets are available already from the ddb
    if isfield(sample_data{k}.meta, 'deployment')
        if isfield(sample_data{k}.meta.deployment, 'StartOffset')
            startOffsets(k) = sample_data{k}.meta.deployment.StartOffset;
        end
        if isfield(sample_data{k}.meta.deployment, 'EndOffset')
            endOffsets(k) = sample_data{k}.meta.deployment.EndOffset;
        end
        if all(isfield(sample_data{k}.meta.deployment, {'TimeDriftInstrument', 'TimeDriftGPS'}))
            if ~isempty(sample_data{k}.meta.deployment.TimeDriftInstrument) && ~isempty(sample_data{k}.meta.deployment.TimeDriftGPS)
                endOffsets(k) = (sample_data{k}.meta.deployment.TimeDriftInstrument - sample_data{k}.meta.deployment.TimeDriftGPS)*3600*24;
            end
        end
    end
end

if ~auto
    f = figure(...
        'Name',        'Time start and end drift calculations in seconds',...
        'Visible',     'off',...
        'MenuBar'  ,   'none',...
        'Resize',      'off',...
        'WindowStyle', 'Modal',...
        'NumberTitle', 'off');
    
    cancelButton  = uicontrol('Style',  'pushbutton', 'String', 'Cancel');
    confirmButton = uicontrol('Style',  'pushbutton', 'String', 'Ok');
    
    setCheckboxes  = [];
    startOffsetFields   = [];
    endOffsetFields = [];
    
    for k = 1:nSample
        
        setCheckboxes(k) = uicontrol(...
            'Style',    'checkbox',...
            'String',   descs{k},...
            'Value',    1, ...
            'UserData', k);
        
        startOffsetFields(k) = uicontrol(...
            'Style',    'edit',...
            'UserData', k, ...
            'String',   num2str(startOffsets(k)));
        
        endOffsetFields(k) = uicontrol(...
            'Style',    'edit',...
            'UserData', k, ...
            'String',   num2str(endOffsets(k)));
        
    end
    
    % set all widgets to normalized for positioning
    set(f,              'Units', 'normalized');
    set(cancelButton,   'Units', 'normalized');
    set(confirmButton,  'Units', 'normalized');
    set(setCheckboxes,  'Units', 'normalized');
    set(startOffsetFields,   'Units', 'normalized');
    set(endOffsetFields,   'Units', 'normalized');
    
    set(f,             'Position', [0.2 0.35 0.6 0.0222 * (nSample + 1)]); % need to include 1 extra space for the row of buttons
    
    rowHeight = 1 / (nSample + 1);
    
    set(cancelButton,  'Position', [0.0 0.0  0.5 rowHeight]);
    set(confirmButton, 'Position', [0.5 0.0  0.5 rowHeight]);
    
    for k = 1:nSample
        
        rowStart = 1.0 - k * rowHeight;
        
        set(setCheckboxes (k), 'Position', [0.0 rowStart 0.6 rowHeight]);
        set(startOffsetFields  (k), 'Position', [0.6 rowStart 0.2 rowHeight]);
        set(endOffsetFields  (k), 'Position', [0.8 rowStart 0.2 rowHeight]);
    end
    
    % set back to pixels
    set(f,              'Units', 'normalized');
    set(cancelButton,   'Units', 'normalized');
    set(confirmButton,  'Units', 'normalized');
    set(setCheckboxes,  'Units', 'normalized');
    set(startOffsetFields,   'Units', 'normalized');
    set(endOffsetFields,   'Units', 'normalized');
    
    % set widget callbacks
    set(f,             'CloseRequestFcn',   @cancelCallback);
    set(f,             'WindowKeyPressFcn', @keyPressCallback);
    set(setCheckboxes, 'Callback',          @checkboxCallback);
    set(startOffsetFields,  'Callback',     @startoffsetFieldCallback);
    set(endOffsetFields,  'Callback',       @endoffsetFieldCallback);
    set(cancelButton,  'Callback',          @cancelCallback);
    set(confirmButton, 'Callback',          @confirmCallback);
    
    set(f, 'Visible', 'on');
    
    uiwait(f);
end

% calculate the drift and apply to the selected datasets
for k = 1:nSample
    
    % this set has been deselected
    if ~sets(k), continue; end
    
    %check for zero values in both fields
    if startOffsets(k) == 0 && endOffsets(k) == 0
        continue;
    end
    
    % look time through dimensions
    type = 'dimensions';
    timeIdx = getVar(sample_data{k}.(type), 'TIME');
    
    if timeIdx == 0
        % look time through variables
        type = 'variables';
        timeIdx = getVar(sample_data{k}.(type), 'TIME');
    end
    
    % no time dimension nor variable in this dataset
    if timeIdx == 0, continue; end
    
    timeDriftComment = ['timeDriftPP: TIME values and time_coverage_end global attributes have been have been '...
        'linearly adjusted for an offset of: ' num2str(startOffsets(k)) ' seconds and a drift of: ' num2str(endOffsets(k) - startOffsets(k)) ' seconds ' ...
        'across the deployment.'];
    
    % apply the drift correction
    newtime = timedrift_corr(sample_data{k}.(type){timeIdx}.data, ...
        startOffsets(k),endOffsets(k));
    sample_data{k}.(type){timeIdx}.data = newtime;
    
    % and to the time coverage atttributes
    sample_data{k}.time_coverage_start = newtime(1);
    sample_data{k}.time_coverage_end = ...
        newtime(end);
    
    
    comment = sample_data{k}.(type){timeIdx}.comment;
    if isempty(comment)
        sample_data{k}.(type){timeIdx}.comment = timeDriftComment;
    else
        sample_data{k}.(type){timeIdx}.comment = [comment ' ' timeDriftComment];
    end
    
    history = sample_data{k}.history;
    if isempty(history)
        sample_data{k}.history = sprintf('%s - %s', datestr(now_utc, ...
            readProperty('exportNetCDF.dateFormat')), timeDriftComment);
    else
        sample_data{k}.history = sprintf('%s\n%s - %s', history, ...
            datestr(now_utc, readProperty('exportNetCDF.dateFormat')), timeDriftComment);
    end
end

    function keyPressCallback(source,ev)
        %KEYPRESSCALLBACK If the user pushes escape/return while the dialog has
        % focus, the dialog is cancelled/confirmed. This is done by delegating
        % to the cancelCallback/confirmCallback functions.
        %
        if     strcmp(ev.Key, 'escape'), cancelCallback( source,ev);
        elseif strcmp(ev.Key, 'return'), confirmCallback(source,ev);
        end
    end

    function cancelCallback(source,ev)
        %CANCELCALLBACK Cancel button callback. Discards user input and closes the
        % dialog .
        %
        sets(:)    = 0;
        startOffsets(:) = 0;
        delete(f);
    end

    function confirmCallback(source,ev)
        %CONFIRMCALLBACK. Confirm button callback. Closes the dialog.
        %
        delete(f);
    end

    function checkboxCallback(source, ev)
        %CHECKBOXCALLBACK Called when a checkbox selection is changed.
        % Enables/disables the offset text field.
        %
        idx = get(source, 'UserData');
        val = get(source, 'Value');
        
        sets(idx) = val;
        
        if val, val = 'on';
        else    val = 'off';
        end
        
        set([startOffsetFields(idx), endOffsetFields(idx)], 'Visible', val);
        
    end

    function startoffsetFieldCallback(source, ev)
        %OFFSETFIELDCALLBACK Called when the user edits one of the offset fields.
        % Verifies that the text entered is a number.
        %
        
        val = get(source, 'String');
        idx = get(source, 'UserData');
        
        val = str2double(val);
        
        % reset the offset value on non-numerical
        % input, otherwise save the new value
        if isnan(val), set(source, 'String', num2str(startOffsets(idx)));
        else           startOffsets(idx) = val;
            
        end
    end

    function endoffsetFieldCallback(source, ev)
        %OFFSETFIELDCALLBACK Called when the user edits one of the offset fields.
        % Verifies that the text entered is a number.
        %
        
        val = get(source, 'String');
        idx = get(source, 'UserData');
        
        val = str2double(val);
        
        % reset the offset value on non-numerical
        % input, otherwise save the new value
        if isnan(val), set(source, 'String', num2str(endOffsets(idx)));
        else           endOffsets(idx) = val;
            
        end
    end

    function newtime = timedrift_corr(time,offset_s,offset_e)
        % remove linear drift of time (in days) from any instrument.
        % the drift is calculated using the start offset (offset_s in seconds) and the
        % end offset (offset_e in seconds).
        
        % calculate the offset times in days:
        offset_days_e = offset_e/60/60/24;
        offset_days_s = offset_s/60/60/24;
        
        if offset_e == offset_s % then just remove the start time
            newtime = time - offset_days_s;
        else
            % make an array of time corrections using the offsets:
            tarray = (offset_days_s:(offset_days_e-offset_days_s)/(length(time)-1):offset_days_e)';
            newtime = time - tarray;
        end
    end
end
