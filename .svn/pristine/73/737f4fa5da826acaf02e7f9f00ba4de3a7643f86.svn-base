function [out,indx] = sort_nd(in)
% function [out {,indx}] = sort_nd(in)
%
%   Sort a vector, removing duplicates. Optionally provide an index.
%
% AUTHOR:   J R Dunn, CSIRO    18/3/96

j = length(in);

if j<2
  out = in;
  indx = 1;
else
  if nargout==1
    out = sort(in);
    remov = find(out(1:j-1)==out(2:j));
    out(remov) = [];
  else
    [out,indx] = sort(in);
    remov = find(out(1:j-1)==out(2:j));
    out(remov) = [];
    indx(remov) = [];
  end
end
