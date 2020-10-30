function datetime = simrad_date_string(filename)
% simrad_date_string extracts a string in the format D20110921-T131520 from
% the filename of a Simrad sounder file.

    [~, name, ~] = fileparts(filename);

    dt = regexp(name, '(\d{8})-T?(\d{6})','tokens');

    if isempty(dt)
        fprintf('Date and time not found in %s\n', name);
        error 'Data file name does not include -D and -T date/time'
    end
        
    datetime = [ 'D' dt{1}{1} '-T' dt{1}{2} ];
end
