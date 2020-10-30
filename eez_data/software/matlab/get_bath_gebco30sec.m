% GET_BATH_GEBCO30SEC  Get GEBCO08 30sec bathymetry at full resolution
%
%  **** See www.marine.csiro.au/eez_data/doc/bathy.html 
%
% INPUTS
%  rgn  - region [w e s n]
%
% OUTPUTS
%  deps   depths (m), +ve upwards, NaN where no value
%  x,y    location vectors
%
% Jeff Dunn   6/6/2012
%
% USAGE: [deps,x,y] = get_bath_gebco30sec(rgn);

function [deps,x,y] = get_bath_gebco30sec(rgn)

fnm = '/home/eez_data/bath/gebco08_30sec.nc';


lo = getnc(fnm,'longitude');
la = getnc(fnm,'latitude');

ix = find(lo>=rgn(1) & lo<=rgn(2));
iy = find(la>=rgn(3) & la<=rgn(4));

x = lo(ix);
y = la(iy);
deps = getnc(fnm,'height',[iy(1) ix(1)],[iy(end) ix(end)]);

%---------------------------------------------------------------------------
