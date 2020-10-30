% COORD2GRD   Translate coordinates onto a grid by scaling, rotation and/or
%             translation
%
% INPUT  lo,la  Matrices of X and Y
%        x0,y0  Rotation and/or translation origin
%        xinc,yinc  Grid spacing in X and Y directions
%        rot    Rotation angle (+ve counter-clock) in radians.
%
% OUTPUT:  x,y  Transformed coordinates (indexed from [1,1] at the origin)
%
%   Author: J Dunn  CSIRO DMR  29/5/98
%
% USAGE: [x,y] = coord2grd(lo,la,x0,y0,xinc,yinc,rot);

function [x,y] = coord2grd(lo,la,x0,y0,xinc,yinc,rot)

if nargin<4
  error('USAGE: [x,y] = coord2grd(lo,la,x0,y0,xinc,yinc,rot)');
end

if nargin<6
  xinc = 1;
  yinc = 1;
end
if nargin<7 | rot==0
  rot = [];
end

if ~isempty(rot)
  c = cos(rot);
  s = sin(rot);
  x = (lo-x0)*c - (la-y0)*s;
  y = (lo-x0)*s + (la-y0)*c;
else
  x = lo-x0;
  y = la-y0;  
end

x = 1 + x/xinc;
y = 1 + y/yinc;

%-------------------------------------------------------------------
