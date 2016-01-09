function [d, a] = calculate_acceleration(pos1, pos2, m1, geocentric)

%% calculate_acceleration Calculate the distance and acceleration between
% two celestial bodies
% [d, a] = calculate_acceleration(pos1, pos2, m1, geocentric)
% calculates the distance between two celestial bodies, bodies, body 1 and 
% body 2, and the magnitude of the gravitational acceleration on body 2 due
% to body 1 
%
% Given: 
% - the position of body 1, pos1 in km
% - the position of body 2 in km
% - mass of body 1, m1 in kg
% - a flag to indicate if the coordinates of body 1 are geocentric
%
% this model returns:
% - the distance between the objects, d (in m)
% - the gravitional acceleration of body 1 due to body 2
%
% Both pos1 and pos2 must be given as Nx3 matrices, where N is the number
% of points to evaluate. The returned values a and d will be Nx1 vectors.

G = 6.67408e-11;                      %Gravitional constant [N m^2 kg^(-2)]

if (geocentric)
    d = sqrt(sum((pos1').^2))*1000;
else
    d = sqrt(sum((pos1' - pos2').^2))*1000;
end

a = m1*G./(d.^2);

