function [T, N] = calculate_undulation(r, az, el, u, coeffs_zonal, R, MU)

%% calculate_undulation Calculate the undulation of a geoid
% [T,N] = calculate_undulation(r, az, el, u, coeffs_zonal, R, MU)
% calculates the undulation of a geoid using Bruns formula.
%
% Given: 
% - a radius, r
% - an azimuth angle, az
% - an elevation angle, el
% - a potential, u
% - a set of even zonal coefficients, coefs_zonal
% - geoid equitorial radius, R
% - geoid gravitational parameter, MU
%
% this model returns:
% - the disturbing potential, T
% - the geoid undulation, N
%
% The radius must be given in meters and the geocentric azimuth and 
% elevation angles in radians. 
% 
% The zonal coefficients must be provided in a Nx6 matrix, with the the
% columns giving in order, n, m, Cnm, Snm. The last 2 columns provide
% standard deviations and in fact need not be included.
% 
% R must be given in meters and MU in km^3 s^-2

v = geopotential(r, az, el, coeffs_zonal);

T = u - v;

N = R^2*T/MU;
