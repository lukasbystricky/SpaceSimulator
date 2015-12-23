function [r, v] = kepler_to_pv(parameters)

%%inputs:
% parameters(1) : semi-major axis, a                     [DU]
% parameters(2) : eccentricity, e                        []
% parameters(3) : inclination, i                         [degrees]
% parameters(4) : argument of periapsis, omega           [degrees]
% parameters(5) : longitude of ascending node, Omega     [degrees]
% parameters(6) : true anomaly, nu                       [degrees]

a = parameters(1);
e = parameters(2);
i = parameters(3)*pi/180;
omega = parameters(4)*pi/180;
Omega = parameters(5)*pi/180;
nu = parameters(6)*pi/180;

% r0 = a*(1 - e^2)/(1 + e*cos(nu));
% 
% q = [r0*cos(nu); r0*sin(nu); 0];
% q_dot = sqrt(a/1 - e^2)*[-sin(nu); e + cos(nu); 0];

p = a*(1 - e^2);
q = p/(1 + e*cos(nu))*[cos(nu); sin(nu); 0];
q_dot = 1/sqrt(p)*[-sin(nu); e + cos(nu); 0];

R1 = [cos(-Omega),  sin(-Omega), 0;
      -sin(-Omega), cos(-Omega), 0;
      0,            0,           1];

R2 = [1, 0 ,       0;
      0, cos(-i),  sin(-i);
      0, -sin(-i), cos(-i)];
  
R3 = [cos(-omega),  sin(-omega), 0;
      -sin(-omega), cos(-omega), 0;
      0,            0,           1]; 
  
r = R1*R2*R3*q;
v = R1*R2*R3*q_dot;