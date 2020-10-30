% GRD2COORD   Translate from a grid to lon,lat coordinates
%
% INPUT  ix,iy    X and Y grid coords
%        lo0,la0  Rotation and/or translation origin 
%        xinc,yinc  Grid spacing in X and Y directions
%        rot      Rotation angle (+ve counter-clock) in radians.
%
% OUTPUT:  lo,la  Transformed coordinates
%
%          Author: J Dunn  CSIRO DMR  29/5/98
%
% NOTE: Assumes the origin is at coord [1,1]  (not [0,0]) 
%
% NOTE: Differs from COORD2GRD in that we first rotate, then scale and 
%       translate. All args are same as for COORD2GRD - we invert them.
%
%  Example: to create coords for a grid -
%    [xx,yy] = meshgrid(1:20,1:20);
%    [Xg,Yg] = grd2coord(xx,yy,lo0,la0,xinc,yinc,rot); 
%
% USAGE: [lo,la] = grd2coord(ix,iy,lo0,la0,xinc,yinc,rot);

function [lo,la] = grd2coord(ix,iy,lo0,la0,xinc,yinc,rot)

if nargin<4
  error('USAGE: [lo,la] = grd2coord(ix,iy,lo0,la0,xinc,yinc,rot)');
end

if nargin<6
  xinc = 1;
  yinc = 1;
end
if nargin<7 | rot==0
  rot = [];
end
 
if ~isempty(rot)
  c = cos(-rot);
  s = sin(-rot);
  lo = (x-1)*c - (y-1)*s;
  la = (x-1)*s + (y-1)*c;
else
  lo = x-1;
  la = y-1;
end

lo = (lo*xinc)+lo0;
la = (la*yinc)+la0;  


%-------------------------------------------------------------------
