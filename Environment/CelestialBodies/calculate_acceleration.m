function [d, a] = calculate_acceleration(pos1, pos2, m1)

G = 6.67408e-11;
d = sqrt(sum((pos1' - pos2').^2))*1000;
a = m1*G./(d.^2);

