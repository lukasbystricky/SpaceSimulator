function [u, g] = geopotential(r, az, el, coeffs)

%% geopotential Spherical harmonic expansion to calculate geopotential and
% accleration
% [u, g] = geopotential(r, az, el, coeffs) implements a spherical harmonic
% expansion to calculate the gravitational potential of a geoid. The
% acceleration is also calculated.
%
% Given: 
% - a radius,r 
% - an azimuth angle, az
% - an elevation angle, el
% - a set of harmonic coefficients, coeffs
% this model returns:
% - the gravitational potential, u
% - the gravitational accleration, g
%
% The radius is given in meters and the geocentric elevation and azimuth 
% angles are given in radians. The coefficients must be fully normalized
% and given in a Nx6 matrix, with the first 4 columns being in order: n, m,
% Cnm, Snm. 
%
% The accleration returned is a 3x1 vector. This vector is first calculated
% as the gradient of the potential in spherical coordinates. It is then
% transformed to cartesian coordinates, relative to the prime meridian.
% This method of calculating the acceleration results in singularities at
% the poles. 

MU = 3.986004415e14;      % Gravitational parameter of Earth    [km^3 s^-2]
A = 6.3781363e6;          % Mean equitorial radius              [m]

u    = MU/r;
g    = zeros(3,1);
g(1) = -MU/r^2;

[nrows, ncols] = size(coeffs);

if (ncols > 1)
    N = max(max(coeffs(:,1:2)));   
    
    [P, Pd] = associated_legendre(N + 1, sin(el), 'egm96');
    
    for i = 1:nrows
        n = coeffs(i,1);
        m = coeffs(i,2);
        C = coeffs(i,3);
        S = coeffs(i,4);
        
        Pnm = P(m+1,n+1);
        Pnmd = Pd(m+1,n+1);
        
        u = u + A^n*MU*(Pnm*(C*cos(m*az) + S*sin(m*az)))/(r^(n+1));
        
        g(1) = g(1) - A^n*(n+1)*MU*Pnm*(C*cos(m*az) + S*sin(m*az))/(r^(n+2));
        g(2) = g(2) + A^n*MU*Pnmd*cos(el)*(C*cos(m*az) + S*sin(m*az))/(r^(n+2));
        g(3) = g(3) + A^n*MU*m*Pnm*(S*cos(m*az) - C*sin(m*az))/(r^(n+2)*sin(el));
    end
end
rot = [sin(el)*cos(az), cos(el)*cos(az), -sin(az);
       sin(el)*sin(az), cos(el)*sin(az), cos(az);
       cos(el),         -sin(el),        0];

g = rot*g;

