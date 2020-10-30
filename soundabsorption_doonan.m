function sa = soundabsorption_doonan(depth, temp, sal, soundSpeed, frequency)
% Calculate sound absorption in dB/m
%
%  Algorithm based on I.J. Doonan, R. F. Coombs and S McClatchie
%  "Absorption of sound in seawater in relation to estimation of deep
%  water fish biomass" ICES Journal of Marine Science 60: 1047-1055,2003
%  With corrections from the authors.

depth = depth / 1000; % depth in km

cold = temp < 20;
a3 = zeros(size(temp));
            
% MgSO4 relaxation component
a2 = 22.19 .* sal .* (1 + 0.017 .* temp) ./ soundSpeed;
f2 = 1.8 .* 10 .^ ( 7 - 1518 ./ ( temp + 273.1));
p2 = exp(-0.176 * depth);

% Freshwater absorption component
a3( cold) = 4.937e-4 - temp( cold) .* (2.590e-5 - temp( cold) .* (9.11e-7 - temp( cold) .* 1.5e-8));
a3(~cold) = 3.964e-4 - temp(~cold) .* (1.146e-5 - temp(~cold) .* (1.45e-7 - temp(~cold) .* 6.5e-10));
p3 = 1 - depth .* (0.0383 - depth .* 4.9e-4);

freq2 = frequency ^ 2;

p2 = repmat(p2',size(temp,1),1);
p3 = repmat(p3',size(temp,1),1);

mgso4 = a2 .* p2 .* f2 ./ (freq2 + f2 .^ 2);
water = a3 .* p3;

sa = freq2 * (mgso4 + water) / 1000;

