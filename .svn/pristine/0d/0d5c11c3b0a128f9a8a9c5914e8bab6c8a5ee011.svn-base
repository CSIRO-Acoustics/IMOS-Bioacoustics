% GETWMO  Given a range or a polygon, get a list of intersected WMOs
% 
% INPUT: range - either [w e s n]
%                or     [x1 y1; x2 y2; x3 y3; ... xn yn]
%
% OUTPUT: vector of WMO square numbers
%
% Jeff Dunn CSIRO Marine Research  Dec 98
%
% USAGE: wmosq = getwmo(range);

function wmosq = getwmo(range)

if min(size(range))==1
  lo1 = range(1);
  lo2 = range(2);
  la1 = range(3);
  la2 = range(4);
else
  lo1 = min(range(:,1));
  lo2 = max(range(:,1));
  la1 = min(range(:,2));
  la2 = max(range(:,2));
end

wmosq = [];

% Set limits in a fraction so do not choose WMO squares just be because
% specified range is along their edge.

ylim = [la1+.0001  la2-.0001];
xlim = [lo1+.0001  lo2-.0001];
ylim = (ceil(ylim/10)*10)-5;
xlim = (ceil(xlim/10)*10)-5;
latgrd = ylim(1):10:ylim(2);
longrd = xlim(1):10:xlim(2);

% If a square range, all squares are in region of interest, but if a polygon
% range, then only use squares which intersect the range.

if min(size(range))==1
  for ii = latgrd
    for jj = longrd
      wmosq =  [wmosq wmo(jj,ii)];
    end
  end
else
  for ii = latgrd
    for jj = longrd
      edge = isintpl(range(:,1),range(:,2),jj+[-5 5 5 -5],ii+[-5 -5 5 5]);
      inner = isinpoly(jj,ii,range(:,1),range(:,2));
      if edge | inner
	wmosq =  [wmosq wmo(jj,ii)];
      end
    end
  end
end

%--------- End of getwmo ----------------------------------------------
