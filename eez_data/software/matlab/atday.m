% ATDAY:  Evaluate a mean field and temporal harmonics at a given day-of-year
%
% INPUT: doy - scalar or vector day-of-year
%        mn  - mean, scalar or same-size as doy, or any shape if doy is scalar
%        an  - complex annual coeffs, same size as mn
%        sa  - [optional] complex semi-annual coeffs, same size as mn
%        dyr - [optional] year-length. Default is 366, but may want to use 365
%              or other to match model output.
%        cap - [optional] minimum value. Best fit seasonal harmonics may sometimes 
%              create -ve values for some part of the year (which for most
%              properties is not physically reasonable.) Eg for nutrients
%              could use 0, for MLD might use 3?, for T use -3.5?
%
% OUTPUT: val - value at day-of-year
%
% Jeff Dunn CSIRO Marine Research   Last mod: 21/11/07
%
% USAGE:  val = atday(doy,mn,an,sa,dyr,cap);

function val = atday(doy,mn,an,sa,dyr,cap)

if nargin==0
  disp('val = atday(doy,mn,an,sa)');
  return
end

if nargin<5 | isempty(dyr)
   dyr = 366;
end
doy = pi*doy/dyr;

if nargin<6
   cap = [];
end

ii = find(isnan(an));
if ~isempty(ii)
  an(ii) = zeros(size(ii));
end

val = mn + real(an.*exp(-i*2*doy));

if nargin>=4 & ~isempty(sa)
  ii = find(isnan(sa));
  if ~isempty(ii)
    sa(ii) = zeros(size(ii));
  end
  val = val + real(sa.*exp(-i*4*doy));
end

if ~isempty(cap)
   val(val<cap) = cap;
end
