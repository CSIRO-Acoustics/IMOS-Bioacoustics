function data = print_sv(ncfile, imagefile, channel, resolution, dotspercell, fontsize)
% PRINT_SV
% function to write a thumbnail image from IMOS-SOOP-BA Sv data 
% usage: data = print_sv(ncfile, imagefile, resolution, dotspercell, fontsize)
%
%   ncfile - IMOS-SOOP-BA NetCDF file name. default ask user - uigetfile
%   imagefile - Image file to create. default: [ncfile '.png']
%   resolution - dots per inch. default: 300
%   dotspercell - dots per cell of data. default:1
%   fontsize - size of font. default: 1500 * dotspercell / resolution
%
% If ncfile is not provided the user will be asked to select it.
% If imagefile is not provided '.png' will be added to the ncfile name.
%
% Author:   Gordon Keith <gordon.keith@csiro.au>
% Version:  1.0
% Date:     2011-09-19
%

    if nargin == 0 || isempty(ncfile)
        [filename pathname] = uigetfile(...
            {'*.nc' '*.nc NetCDF file'; '*' 'All files'}, ...
            'IMOS-SOOP-BA NetCDF file');
        if isequal(filename,0)
            error('NetCDF filename required');
        end
        ncfile = fullfile(pathname,filename);
    end
    
    if nargin < 2
        imagefile = [ ncfile '.png' ];
    end

    if nargin < 3
        channel = '38';
    end

    if nargin < 4
        resolution = 300;
    end

    if nargin < 5
        dotspercell = 1;
    end

    if nargin < 6
        fontsize = 1500 * dotspercell / resolution;
    end

    if isnumeric(channel)
        channel = num2str(channel);
    end
    
% open the netcdf file
    if isempty(ls(ncfile))
        fprintf('\n-----------------------------------------------------\n');
        fprintf('\nnetcdf file %s \ncannot be found, check location\n',ncfile);
        fprintf('\n-----------------------------------------------------\n');
    else
        ncid = netcdf.open(ncfile, 'NC_NOWRITE');            
        data.file = ncfile;
        [~, filename, ext] = fileparts(ncfile);
        
        try 
            start_location = netcdf.getAtt(ncid, ...
                netcdf.getConstant('NC_GLOBAL'), 'transit_start_locality');
        catch e
            start_location = '';
        end
        try 
            end_location = netcdf.getAtt(ncid, ...
                netcdf.getConstant('NC_GLOBAL'), 'transit_end_locality');
        catch e
            end_location = '';
        end
            
        % read data
        try
            depthid =  netcdf.inqVarID(ncid, 'DEPTH');
        catch e
            depthid =  netcdf.inqVarID(ncid, 'RANGE');
        end
        depth = netcdf.getVar(ncid, depthid);

        timeid =  netcdf.inqVarID(ncid, 'TIME');                    
        time50 = netcdf.getVar(ncid, timeid);
        time = time50 + datenum('1950 01 01');
        
        svid = netcdf.inqVarID(ncid, ['Sv_' channel]);    
        data.sv = netcdf.getVar(ncid, svid);
        data.Sv = 10*log10(data.sv);
                               
        qcid = netcdf.inqVarID(ncid, ['Sv_' channel '_quality_control']);    
        data.qc = netcdf.getVar(ncid, qcid);

        % Ignore bad data
        data.Sv(data.qc>2)=NaN;     
         
        % viz the sv data
        figure
        
        % size figure to hold full data set
        dwidth = size(data.sv,2);
        wdth = dwidth * dotspercell / resolution;
        hght = size(data.sv,1) * dotspercell / resolution;
        left = 0.5 + 50 / resolution;
        bot  = 4 * fontsize / 72;
        width = wdth + 2 * left;
        height = hght + 2 * bot;
        set(gcf, 'Units', 'inches');
        set(gcf, 'Position', [0, 0.5, width, height ] );
        
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [width height]);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [0 0 width height]);
        
        set(gca, 'FontSize', fontsize);
        
        min_sv = -84; range = 36; max_sv = min_sv + range;
        imagesc(data.Sv, [min_sv max_sv])
        colormap(EK500colourmap())
        
        xlabel('Time UTC')
        ylabel('Depth (m)')
        title({ [ filename ext ] ;  'Sv mean (dB re 1 m-1)' }, ...
            'Interpreter','none')
        colorbar('FontSize', fontsize, 'YTick', -80:10:-50);
        
        text(0, -15, start_location, 'FontSize', fontsize);
        text(dwidth, -15, end_location, 'HorizontalAlignment', 'right', 'FontSize', fontsize);

        set(gca, 'Units', 'inches');
        set(gca, 'Position', [left, bot, wdth, hght]);

        ticks(gcf)
        drawnow;
        
        print(gcf, '-dpng', imagefile, ['-r' num2str(resolution)]);
        
        close;
    end


function [EK500cmap]=EK500colourmap()

EK500cmap = [255   255   255  % white
           159   159   159    % light grey
           95    95    95     % grey
           0     0   255      % dark blue
           0     0   127      % blue
           0   191     0      % green
           0   127     0      % dark green
           255   255     0    % yellow
           255   127     0    % orange
           255     0   191    % pink
           255     0     0    % red
           166    83    60    % light brown
           120    60    40]./255;  % dark brown
end

function ticks(figure,~)
    axes=get(figure,'CurrentAxes');
    
    % depth ticks
    ytick=get(axes,'YTick');
    set(axes,'YTickLabel',depth(floor(ytick)));
    
    % time ticks
    % xtick=get(axes,'XTick');
    % set(axes,'XTickLabel',datestr(time(floor(xtick))));
    xlim=get(axes,'XLim');
    start=time(ceil(xlim(1)));
    finish=time(floor(xlim(2)));
    
    len=finish-start;
    if len > 2
        format=29;      % 'yyyy-mm-dd'
        tock=1;         % day
    elseif len > .5
        format = 'yyyy-mm-dd HH:MM';
        tock=4;         % 6 hr
    elseif len > .1
        format = 15;    % 'HH:MM'
        tock=24;        % hour
    else
        format = 13;    % 'HH:MM:SS'
        tock=96;        % quarter hour
    end
    
    xtock=(ceil(start*tock):1:finish*tock)/tock;
    xtick(length(xtock))=0;
    for i=1:length(xtock)
        xtick(i)=find(time>xtock(i),1);
    end
    
    set(axes,'XTick',xtick);
    set(axes,'XTickLabel',datestr(time(xtick), format));
end
end
