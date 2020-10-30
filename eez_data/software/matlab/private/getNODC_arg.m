function [] = getNODC_arg( arg );

% Used by getNODC_var to decode input arguments. There must be a better way!
%
%  Copyright (C) J R Dunn, CSIRO, 
%  $Revision: 1.6 $    Last revision $Date: 1997/02/12 03:34:28 $

global getN_varpr; global getN_varnm; global getN_deps; global getN_getd;
global getN_scr; global getN_nvar;


if getN_getd 
  if isstr(arg)
    error(['Expecting a depth range instead of argument ' arg]);
  else
    getN_getd = 0;
    getN_deps = [getN_deps; arg];
  end
else
  if ~isstr(arg)
    error(['Expecting a variable type string in argument list.']);
  else
    getN_getd = 1;
    getN_nvar = getN_nvar+1;

    if strcmp(arg,'t')
      getN_varpr = 'ts';
      getN_varnm = str2mat(getN_varnm,'t');
    elseif strcmp(arg,'ts')
      getN_varpr = 'ts';
      getN_varnm = str2mat(getN_varnm,'t');
      getN_scr(getN_nvar) = 1;
    elseif strcmp(arg,'ft')
      getN_varpr = 'ts';
      getN_varnm = str2mat(getN_varnm,'t_flag');

    elseif strcmp(arg,'s')
      getN_varpr = 'ts';
      getN_varnm = str2mat(getN_varnm,'s');
    elseif strcmp(arg,'ss')
      getN_varpr = 'ts';
      getN_varnm = str2mat(getN_varnm,'s');
      getN_scr(getN_nvar) = 1;
    elseif strcmp(arg,'fs')
      getN_varpr = 'ts';
      getN_varnm = str2mat(getN_varnm,'s_flag');

    elseif strcmp(arg,'o')
      getN_varpr = 'o2';
      getN_varnm = str2mat(getN_varnm,'o2');
    elseif strcmp(arg,'os')
      getN_varpr = 'o2';
      getN_varnm = str2mat(getN_varnm,'o2');
      getN_scr(getN_nvar) = 1;
    elseif strcmp(arg,'fo')
      getN_varpr = 'o2';
      getN_varnm = str2mat(getN_varnm,'o2_flag');

    elseif strcmp(arg,'n')
      getN_varpr = 'no3';
      getN_varnm = str2mat(getN_varnm,'no3');
    elseif strcmp(arg,'ns')
      getN_varpr = 'no3';
      getN_varnm = str2mat(getN_varnm,'no3');
      getN_scr(getN_nvar) = 1;
    elseif strcmp(arg,'fn')
      getN_varpr = 'no3';
      getN_varnm = str2mat(getN_varnm,'no3_flag');

    elseif strcmp(arg,'p')
      getN_varpr = 'po4';
      getN_varnm = str2mat(getN_varnm,'po4');
    elseif strcmp(arg,'ps')
      getN_varpr = 'po4';
      getN_varnm = str2mat(getN_varnm,'po4');
      getN_scr(getN_nvar) = 1;
    elseif strcmp(arg,'fp')
      getN_varpr = 'po4';
      getN_varnm = str2mat(getN_varnm,'po4_flag');

    elseif strcmp(arg,'i')
      getN_varpr = 'si';
      getN_varnm = str2mat(getN_varnm,'si');
    elseif strcmp(arg,'is')
      getN_varpr = 'si';
      getN_varnm = str2mat(getN_varnm,'si');
      getN_scr(getN_nvar) = 1;
    elseif strcmp(arg,'fi')
      getN_varpr = 'si';
      getN_varnm = str2mat(getN_varnm,'si_flag');

    elseif strcmp(arg,'g')
      getN_varnm = str2mat(getN_varnm,'neut_density');

    elseif strcmp(arg,'cp')
      getN_varnm = str2mat(getN_varnm,'csiro_profile_no');
      getN_getd = 0;
      getN_deps = [getN_deps; 0 0];
    elseif strcmp(arg,'cf')
      getN_varnm = str2mat(getN_varnm,'csiro_flag');
      getN_getd = 0;
      getN_deps = [getN_deps; 0 0];
    elseif strcmp(arg,'cr')
      getN_varnm = str2mat(getN_varnm,'cruise_no');
      getN_getd = 0;
      getN_deps = [getN_deps; 0 0];
    elseif strcmp(arg,'dp')
      getN_varnm = str2mat(getN_varnm,'number_of_depths');
      getN_getd = 0;
      getN_deps = [getN_deps; 0 0];
    elseif strcmp(arg,'co')
      getN_varnm = str2mat(getN_varnm,'country_code');
      getN_getd = 0;
      % Country code is a 2 element variable:
      getN_deps = [getN_deps; 1 2];

    else
      disp(['### Do NOT understand code ' arg]);
      getN_nvar = getN_nvar-1;
    end

  end
end
