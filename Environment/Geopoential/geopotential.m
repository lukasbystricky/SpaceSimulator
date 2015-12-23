function [u, g] = gravitational_acceleration_harmonics(r, az, el, coeffs)

MU = 3.986004415e14;
A = 6.3781363e6;

u    = -MU/r;
g    = zeros(3,1);
g(1) = MU/r^2;

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
        
        u = u - A^n*MU/r*(Pnm*(C*cos(m*az) + S*sin(m*az)))/(r^(n+1));
        
        g(1) = g(1) + A^n*(n+1)*MU*Pnm*(C*cos(m*az) + S*sin(m*az))/(r^(n+2));
        g(2) = g(2) - A^n*MU*Pnmd*cos(el)*(C*cos(m*az) + S*sin(m*az))/(r^(n+2));
        g(3) = g(3) - A^n*MU*m*Pnm*(S*cos(m*az) - C*sin(m*az))/(r^(n+2)*sin(el));
    end
end
rot = [sin(el)*cos(az), cos(el)*cos(az), -sin(az);
       sin(el)*sin(az), cos(el)*sin(az), cos(az);
       cos(el),         -sin(el),        0];

g = rot*g;

