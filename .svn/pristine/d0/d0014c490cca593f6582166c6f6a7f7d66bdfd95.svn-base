function ss = soundspeed_coppens(depth, temp, sal, latitude)
%  Calculate the speed of sound in m/s.
%
% Algorithm based on Kongsberg Simrad EM300 Operator Manual pp376-377,
% based on Coppens, JASA March 1981.
% For depths greater than 5000 m a latitude correction is required, the
% latitude parameter is only required if the depth parameter extends that
% far.

depth = depth / 1000; % depth in km

ss = 1449.05 + temp .* (4.57 - temp .* ( 0.0521 - temp .* 0.00023)) + ...
    ( 1.333 - temp .* (0.0126 - temp .* 0.00009)) .* (sal - 35);

for d = 1 : length(depth)
    if depth(d) <= 0
    elseif depth(d) < 1
        ss(:,d) = ss(:,d) + depth(d) * 16.5;
    elseif depth(d) < 5
        ss(:,d) = ss(:,d)  + depth(d) * ...
            (16.3 + depth(d) * (0.22 - depth(d) * 0.003 * sqrt(max(0,temp(:,d) + 2))));
    else
        ss(:,d) = ss(:,d)  + depth(d) * (1 - 0.0026 .* cos (latitude * 4 * pi / 180)) .* ...
            (16.3 + depth(d) * (0.22 - depth(d) * 0.003 * sqrt(max(0,temp(:,d) + 2))));        
    end
end
