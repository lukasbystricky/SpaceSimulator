function [st_true, st_mean] = calculate_solar_time(gmt, az, d)
% gmt  - hour in day (decimal, Greenwith Mean Time)
% az   - azimuth angle, equivalent to longitude (radians)
% d    - day in year (integer)

st_mean = gmt + az*2*pi/24;

eot = equation_of_time(d);

st_true = st_mean + eot/60;

end

function eot = equation_of_time(d)

B = (360/365)*(d - 81);
eot = 9.87*sind(2*B) - 7.53*cosd(B) - 1.5*sind(B);
end