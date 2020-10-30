function ss = soundspeed_gsw(depth, temp, salt, latitude, longitude)
% calculate sound speed using the Gibbs SeaWater (GSW) oceanographic
% toolbox of TEOS-10 http://www.TEOS-10.org

depth = depth';
height = - depth(ones(size(latitude)),:);
p = gsw_p_from_z(height,latitude);
[SA, ~] = gsw_SA_from_SP(salt,p,longitude,latitude);
CT = gsw_CT_from_t(SA,temp,p);
ss = gsw_sound_speed(SA,CT,p);