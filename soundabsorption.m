function sa = soundabsorption(depth, temp, sal, soundspeed, frequency, pH)
% Calculate sound absorption in dB/m
%
% Algorithm based on Kongsberg Simrad EM300 Operator Manual pp378-379,
% based on R. E. Francois and G.R. Garrison "Sound absorption based on ocean measurements:
% Part II: Boric acid contribution and equation for total absorption,"
% J of Acoust Soc Am 72(6) Dec 1982 p 1886

depth = depth / 1000; % depth in km

cold = temp < 20;
a3 = zeros(size(temp));

a1 = 8.86 * 10 .^ (0.78 * pH - 5) ./ soundspeed;
a2 = 21.44 .* sal .* (1 + 0.025 .* temp) ./ soundspeed;
a3(cold)  = 4.937e-4 - temp( cold) .* (2.590e-5 - temp( cold) .* (9.11e-7 - temp( cold) * 1.5e-8));
a3(~cold) = 3.964e-4 - temp(~cold) .* (1.146e-5 - temp(~cold) .* (1.45e-7 - temp(~cold) * 6.5e-10));

p2 = 1 - depth .* (0.137 - depth * 0.0062);
p3 = 1 - depth .* (0.0383 - depth * 4.9e-4);

f1 = 2.8 * sqrt(sal / 35) .* 10 .^ (4 - 1245 ./ (273 + temp));
f2 = 8.17 * 10 .^ (8 - 1990 ./ (273 + temp)) ./ (1 + 0.0018 * (sal - 35));

freq2 = frequency ^ 2;

mgso4 = zeros(size(temp));
water = zeros(size(temp));

boric = a1 .* f1 ./ (freq2 + f1 .* f1);
for i = 1 : size(temp,1);
    mgso4(i,:) = a2(i,:) .* p2' .* f2(i,:) ./ (freq2 + f2(i,:) .^ 2);
    water(i,:) = a3(i,:) .* p3';
end

sa = freq2 * (boric + mgso4 + water) / 1000;
