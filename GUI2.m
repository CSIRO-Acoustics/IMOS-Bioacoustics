function GUI2(defaults)
% Call GUI with the NetCDF mode.
% This function is just a shortcut to GUI(defaults, 'NetCDF')
%
% The NetCDF mode contains the GUI settings for all processing from reading
% echointegration results.

if nargin < 1
    defaults = '?';
end

GUI(defaults, 'NetCDF');