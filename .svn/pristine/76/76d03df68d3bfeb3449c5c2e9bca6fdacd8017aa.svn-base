function [var] = getNODC_onev(netcdf_file,nv,var,reject,nd,nc);

%  Used by getNODC_var to read one variable, and append it to the data 
%  collected so far.
%
% INPUTS:
%  netcdf_file - name of file
%  nv   - index of this one variable (in the global arrays). 
%  var  - data already accumulated
%  reject - list of indices of casts to be rejected
%  nd   - number of depths in this file
%  nc   -   "    "  casts  "       "
%
%  Copyright (C) J R Dunn, CSIRO, 
%  $Revision: 1.3 $    Last revision $Date: 1997/11/21 02:43:16 $

% Globals:
%  varnm - names of variables in netcdf file
%  deps -  array of [d1 d2] specifying (upper and lower) range of standard 
%          depths to get

global getN_varnm; global getN_deps; global getN_scr;


vntmp = deblank(getN_varnm(nv,:));

if getN_deps(nv,:)==[0 0]
  % This should only be used to extract 1-D variables
  VAR = getcdf(netcdf_file,vntmp);
elseif getN_deps(nv,1) <= nd
  % getcdf crashes if ask for more depths than are in the file - so work out
  % max depth to ask for:
  maxdep = min([nd getN_deps(nv,2)]);
  VAR = getcdf(netcdf_file,vntmp,[-1 getN_deps(nv,1)],[1 maxdep]);

  % We want each row to correspond to a cast, so transpose. However, if only
  % extract one depth, it already comes out as a column-cast vector, so
  % don't change.
  if (maxdep>getN_deps(nv,1));  VAR = VAR'; end

  % If required, screen using individual data flags.
  if getN_scr(nv)
    vntmp = [vntmp '_flag'];
    flags = getcdf(netcdf_file,vntmp,[-1 getN_deps(nv,1)],[1 maxdep]);
    if (maxdep>getN_deps(nv,1));  flags = flags'; end
%#    rejs = find(flags~=0);
    rejs = find(rem(flags,10)~=0);
    VAR(rejs) = NaN*ones(size(rejs));
  end
  
else
  % Case where the specified depth levels are not present in the file,
  % so make up a dummy VAR.
  VAR = NaN*ones(nc,1);
end


VAR(reject,:) = [];


% If necessary, pack the data to the required number of depths.

[m n] = size(VAR);
ndeps = 1+getN_deps(nv,2)-getN_deps(nv,1);
if n < ndeps
  VAR = [VAR NaN*ones(m,ndeps-n)];
end

var = [var; VAR];

% --------------- End of getNODC_onev.m ------------------
