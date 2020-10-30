function meta = read_ecs(ecs_file)
% Read an echoview extended calibration file and write the calibration
% settings to the metadata structure.
%
% Inputs:
%   ecs_file    full path to EchoView Calibration Settings file
%
% Outputs:
%   meta        structure containing metadata
%
% An ecs file contains default settings and active settings. 
% Active settings are in the form "parameter = value" optionally followed
% by a comment ("#" followed by text).
% Default settings are similar to active settings but a preceded by a
% comment character "#".
%
% There are usually multiple sections in a file, each parameter may appear
% in each section.
%
% This function returns the last value of the parameter to appear in the
% file, with active settings taking precedence over default settings
%

active = struct;
defalt = struct;
settings = struct;
cal={};
cals=0;
cal_date = [];
try
    % read cal file
    fid = fopen(ecs_file, 'rt');
    
    line = fgetl(fid);
    while ~isempty(line) && line(1) > 128   % remove UTF code on first line
        line = line(2:end);
    end
    while (ischar(line))
        eq = find(line=='=',1);
        if strncmp('Version', line, 7)
        elseif isempty(line)
        elseif line(1) ~= ' ' && line(1) ~= '#'
            % new calibration section
            if cals 
                % save previous cal
                fields = fieldnames(defalt);
                for field = fields'
                    settings(cals).(field{1}) = defalt.(field{1});
                end
                fields = fieldnames(active);
                for field = fields'
                    settings(cals).(field{1}) = active.(field{1});
                end
            end
            % create new cal
            cals=cals+1;
            cal{cals} = line;           %#ok<AGROW>
            settings(cals).calibration_name = line;
            active = struct;
            defalt = struct;
        elseif eq > 1
            % calibration parameter
            comment = find(line(eq:end)=='#',1);
            if comment
                line=line(1:eq+comment-2);
            end
            
            calval = strtrim(line(eq+1:end));
            if ~isempty(calval)
                numval = str2double(calval);
                if ~isnan(numval)
                    calval = numval;
                end
                
                calparam=strtrim(line(1:eq -1));
                
                % commented parameters are treated as defaults, only used
                % if there are no uncommented versions of the that param
                comment = find(calparam=='#',1);
                if comment
                    calparam=strtrim(calparam(comment + 1:end));
                    if ~isempty(calparam)
                        defalt.(calparam) = calval;
                    end
                else
                    active.(calparam) = calval;
                end
            end
        elseif line(1) == '#' && ~isempty(find(line == '/',1))
            % ecs file date, not sure if its meaningful
            try
                line(line=='#') = '';
                cal_date = datenum(line,'dd/mm/yy HH:MM:SS');
            catch
            end
        end
        
        line = fgetl(fid);
    end    
    
    fclose(fid);
    
    % save last cal
    if cals
        fields = fieldnames(defalt);
        for field = fields'
            settings(cals).(field{1}) = defalt.(field{1});
        end
        fields = fieldnames(active);
        for field = fields'
            settings(cals).(field{1}) = active.(field{1});
        end
    end
    
catch e
    warning('read_ecs:exception', 'Problem reading cal file %s\n%s', ...
        ecs_file, e.message)
end

if cals == 0
    error('No calibrations found in settings file: %s', ecs_file);
end
    
% translate probably should be read from a config file.
% column 1 is ECS parameter, column 2 is metadata field.
% ECS parameter may appear more than once. metdata field may be blank.
% translate = {   
%     'AbsorptionCoefficient'     , 'data_processing_absorption';
%     'EK60SaCorrection'          , 'data_processing_sa_correction';
%     'Ek60TransducerGain'        , 'data_processing_transceiver_gain';
%     'Frequency'                 , 'data_processing_frequency';
%     'MajorAxis3dbBeamAngle'     , 'data_processing_transducer_beam_angle_major';
%     'MajorAxisAngleOffset'      , '';
%     'MajorAxisAngleSensitivity' , '';
%     'MinorAxis3dbBeamAngle'     , 'data_processing_transducer_beam_angle_minor';
%     'MinorAxisAngleOffset'      , '';
%     'MinorAxisAngleSensitivity' , '';
%     'SoundSpeed'                , 'data_processing_soundspeed';
%     'TransmittedPower'          , 'data_processing_transceiver_power';
%     'TransmittedPulseLength'    , 'data_processing_transmit_pulse_length';
%     'TvgRangeCorrection'        , '';
%     'TwoWayBeamAngle'           , 'data_processing_transducer_psi';
%     
%     'AccuracyEstimate'          , 'calibration_accuracy_estimate';
%     'AcquisitionMethod'         , 'calibration_acquisition_method';
%     'CalibrationDate'           , 'calibration_date'
%     'ProcessingMethod'          , 'calibration_processing_method';
%     'ReportURL'                 , 'calibration_report';
%     'TransducerDepth'           , 'data_processing_transducer_depth'
%     };

% edit by Haris to match with updated metadata record (17 November 2017)

translate = {   
    'AbsorptionCoefficient'     , 'data_processing_absorption';
    'EK60SaCorrection'          , 'data_processing_Sacorrection';
    'Ek60TransducerGain'        , 'data_processing_on_axis_gain';
    'Frequency'                 , 'data_processing_frequency';
    'MajorAxis3dbBeamAngle'     , 'data_processing_transducer_beam_angle_major';
    'MajorAxisAngleOffset'      , '';
    'MajorAxisAngleSensitivity' , '';
    'MinorAxis3dbBeamAngle'     , 'data_processing_transducer_beam_angle_minor';
    'MinorAxisAngleOffset'      , '';
    'MinorAxisAngleSensitivity' , '';
    'SoundSpeed'                , 'data_processing_soundspeed';
    'TransmittedPower'          , 'data_processing_transceiver_power';
    'TransmittedPulseLength'    , 'data_processing_transmit_pulse_length';
    'TvgRangeCorrection'        , '';
    'TwoWayBeamAngle'           , 'data_processing_transducer_psi';
    
    'AccuracyEstimate'          , 'calibration_accuracy_estimate';
    'AcquisitionMethod'         , 'calibration_acquisition_method';
    'CalibrationDate'           , 'calibration_date'
    'ProcessingMethod'          , 'calibration_processing_method';
    'ReportURL'                 , 'calibration_report';
    'TransducerDepth'           , 'data_processing_transducer_depth'
    };


% calibration metadata
[~, cfile, cext] = fileparts(ecs_file);
for c = cals:-1:1
    meta(c).calibration_name = cal{c};
    meta(c).calibration_file_name = [cfile cext];
    if ~isempty(cal_date)
        % meta(c).calibration_file_date = cal_date;
    end
end

% translate cal settings to meta data
for i = 1:size(translate,1)
    if ~isempty(translate{i,2})
        for c = cals:-1:1
            if isfield(settings(c), translate{i,1})
                meta(c).(translate{i,2}) = settings(c).(translate{i,1});
            end
        end
    end
end


