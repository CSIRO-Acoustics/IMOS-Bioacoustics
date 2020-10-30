% COMBINE  Combine two lists, removing duplicates, and indexing contributions
%          from both lists.
%
% dtot is union of inputs (with duplicates removed)
% i1 is location of each of l1 in dtot 
% i2 is location of each of l2 in dtot 
% ovl: -1:  No overlap
%       0:  l2 has values not in l1
%       1:  l2 is a subset of l1
%       2:  l2 = l1
%
% Jeff Dunn  30/3/99
%
% USAGE: [ltot,i1,i2,ovlap] = combine(l1,l2);

function [dtot,i1,i2,ovlap] = combine(d1,d2)

n1 = length(d1);
n2 = length(d2);

[dtot,idx] = sort([d1(:); d2(:)]);

dup = [0; (abs(dtot(1:end-1)-dtot(2:end))<.01)];
dtot(find(dup)) = [];

if ~any(dup)
   ovlap = -1; 				% NO overlap between new & old
elseif length(dtot) == n1
   if length(dtot) == n2
      ovlap = 2;
   else
      ovlap = 1;
   end
else
   ovlap = 0;
end

i1 = [];
i2 = [];
jj = 0;
for ii = 1:length(idx)
   if ~dup(ii)
      jj = jj+1;
   end
   if idx(ii)>n1
      i2(idx(ii)-n1) = jj;
   else
      i1(idx(ii)) = jj;
   end      
end

return
%---------------------------------
