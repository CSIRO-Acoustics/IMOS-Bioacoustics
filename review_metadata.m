function sample_data = review_metadata(sample_data)
% Display the IMOS toolbox metadata editing window and all the user to
% update the metadata.
%
    if isempty(which('viewMetadata'))
        addpath(fullfile(fileparts(mfilename('fullpath')), 'IMOS-toolbox', 'GUI'));
        addpath(fullfile(fileparts(mfilename('fullpath')), 'IMOS-toolbox', 'Util'));
        addpath(fullfile(fileparts(mfilename('fullpath')), 'IMOS-toolbox', 'NetCDF'));
        addpath(fullfile(fileparts(mfilename('fullpath')), 'IMOS-toolbox', 'IMOS'));
        addpath(fullfile(fileparts(mfilename('fullpath')), 'IMOS-toolbox', 'Parser'));
    end
    if isempty(which('viewMetadata'))
        addpath(genpath(fileparts(mfilename('fullpath'))));
    end
    
    wd = cd(fileparts(fileparts(which('viewMetadata'))));
    
    fig = figure;    
    viewMetadata(fig, sample_data, @sync, @donothing, 'timeSeries')
    waitfor(fig)
    
    cd(wd);
    
    function sync (sd)
        sample_data = sd;
    end

    function donothing(~,~,~)
        
    end
end