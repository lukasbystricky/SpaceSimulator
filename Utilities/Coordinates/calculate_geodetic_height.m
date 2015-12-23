function h = calculate_geodetic_height(r_geocentric)

R_EQ = 6.378145e6;
R_POL = 6.356785e6;

el = acos(r_geocentric(3)/norm(r_geocentric));

r = sqrt(((R_EQ^2*cos(el))^2 + (R_POL^2*sin(el))^2)/((R_EQ*cos(el))^2 + (R_POL*sin(el))^2));
h = norm(r_geocentric) - r;
