function [T, N] = calculate_undulation(r, az, el, u, coeffs_zonal, R, MU)

v = geopotential(r, az, el, coeffs_zonal);

T = u - v;

N = R^2*T/MU;
