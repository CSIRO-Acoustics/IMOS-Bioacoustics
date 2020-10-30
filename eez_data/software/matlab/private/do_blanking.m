% DO_BLANKING:  Set map values to nan if that point is below ocean bottom
%
%   in:    A [101,161] OceansEEZ map
%   dep:   The standard depth of that map (ie between 1 and 33)
% 
% JRD 20/8/97

function [out] = do_blanking(in,dep,Xg,Yg)

out = repmat(nan,size(in));

load /home/eez_data/atlas/deepest;
[dpx,dpy]=meshgrid(dpx,dpy);
dpg = interp2(dpx,dpy,deepest,Xg,Yg);
ii = find(dpg>=dep);
out(ii) = in(ii);
