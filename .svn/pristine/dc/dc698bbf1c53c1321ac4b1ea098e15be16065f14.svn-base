
%Intergrate to defined depth

setting.Grid_height=5
setting.Grid_distance=500

GUI(setting)

%% Generate a netCDF to defined max_depth 
% Note: that this will effect all frequencies

setting.max_depth=400
GUI(setting)

%% Plot to defined depth

data=viz_sv                         % nav to netCDF file of interest

data4=cut_depth(data,400)           % define depth to plot to
viz_sv(data4,data4.Sv,'range',[])  % plot default Sv vbalue if multi channel then 
viz_sv(data4,data4.Sv,'range',[], 'channel', 120, 'sun') % plot 120kHz channel to 400
viz_sv(data,data.Sv,'range',[], 'channel', 38, 'sun') % plot 38kHz channel to full extent