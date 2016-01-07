function [t_inf, T, c_H, c_He, c_O, c_N2] = dtm94(z, coeffs, P, f, f_mean, k, t, d)

%% dtm94 Drag Temperature Model for high altitude atmospheric density 
% [t_inf, T, cH, c_He, c_O, c_N2] = dtm94(z, coeffs, P, f, f_mean, k, t, d)
% implements the DTM 94 model for high altiude density as desribed by
% Berger et al. in "Improvement of the empirical thermospheric model DTM:
% DTM94 - a comparative review of various temporal variations
% and prospects in space geodesy applications".
%
% Given: 
% - an altiude, z
% - a set of coefficients, coeffs
% - a vector of Legendre polynomials evaluated at sin(latitude), P
% - the daily solar flux, f
% - the average solar flux, f_mean
% - the current geomagnetic index, k
% - the local solar time, t
% - the day number, d
% this model returns:
% - the thermopause temperature, t_inf
% - the temperature, T
% - the concentration of Hydrogen, c_H
% - the concentration of Helium, c_He
% - the concentration of Oxygen, c_O
% - the concentration of molecular Nitrogen, c_N2
%
% The altitude z must be a scalar given in km.
%
% The coefficients must be given in a 36x4 matrix in the order of
% temperature, He, O, N2.
%
% P must be a matrix of size 6x6, such that P(i,j) = the associated
% legendre polynomial of degree i-1 and order j-1 evaluated at
% sin(latitude).
%
% f and f_mean must be measured in 10^-22 W m^-2 Hz^-1.
%
% k must be an integer between 0 and 9.
%
% t must be measured in hours, between 0 and 24.
%
% d is the day number in the year, between 0 and 365.
%
% The temperatures returned are given in degrees Kelvin, while the
% concentrations are given in particles/cm^3.


T120 = 380;              % Temperature at 120 km                  [K]
R = 6.35677e3;           % Equitorial radius of Earth             [km]
T120p = 14.348;          % Relative temperature vertical gradient [K/km]

MU = 398600.4418;        % Gravitational parameter of Earth     [km^3 s^-2]
BOLTZ = 1.38064852;      % Boltzmann constant                   [10^29 km^2 kg s^-2 K-1]

MASS_H = 1.6737236e2;    % Mass of Hydrogen                     [10^29 kg]
MASS_HE = 6.6464764e2;   % Mass of Helium                       [10^29 kg]
MASS_O = 2.6567626e3;    % Mass of Oxygen                       [10^29 kg]
MASS_N2 = 2*2.3258671e3; % Mass of di-Nitrogen                  [10^29 kg]

g = MU/((R+120)^2);

xi = (z - 120)*(R + 120)/(z + R);

[t_inf, cH_0, cHe_0, cO_0, cN2_0] = calculate_t_inf_and_initial_concentrations(...
    coeffs, P, f, f_mean, k, t, d);

sigma = T120p/(t_inf - T120);

gammaH = MASS_H*g/(sigma*BOLTZ*t_inf);
gammaHe = MASS_HE*g/(sigma*BOLTZ*t_inf);
gammaO = MASS_O*g/(sigma*BOLTZ*t_inf);
gammaN2 = MASS_N2*g/(sigma*BOLTZ*t_inf);

alphaH = -0.38;
alphaHe = -0.38;

T = t_inf - (t_inf - T120)*exp(-sigma*xi);

fH = (T120/T)^(1 + alphaH + gammaH)*exp(-sigma*gammaH*xi);
fHe = (T120/T)^(1 + alphaHe + gammaHe)*exp(-sigma*gammaHe*xi);
fO = (T120/T)^(1 + gammaO)*exp(-sigma*gammaO*xi);
fN2 = (T120/T)^(1 + gammaN2)*exp(-sigma*gammaN2*xi);

c_H = cH_0*fH;
c_He = cHe_0*fHe;
c_O = cO_0*fO;
c_N2 = cN2_0*fN2;

end

function [t_inf, cH, cHe, cO, cN2] = calculate_t_inf_and_initial_concentrations(...
    coeffs, P, f, f_mean, k, t, d)

coeffsT = coeffs(:,1);
coeffsH = coeffs(:,2);
coeffsHe = coeffs(:,3);
coeffsO = coeffs(:,4);
coeffsN2 = coeffs(:,5);

GT = harmonic_expansion(coeffsT, P, f, f_mean, k, t, d, false, true);
GH = harmonic_expansion(coeffsH, P, f, f_mean, k, t, d, true, false);
GHe = harmonic_expansion(coeffsHe, P, f, f_mean, k, t, d, true, false);
GO = harmonic_expansion(coeffsO, P, f, f_mean, k, t, d, false, false);
GN2 = harmonic_expansion(coeffsN2, P, f, f_mean, k, t, d, false, false);

t_inf = coeffsT(1)*GT;
cH = coeffsH(1)*exp(GH - 1);
cHe = coeffsHe(1)*exp(GHe - 1);
cO = coeffsO(1)*exp(GO - 1);
cN2 = coeffsN2(1)*exp(GN2 - 1);

end

function G = harmonic_expansion(A, P, f, f_mean, k, t, d, isH, isT)

omega = 2*pi/24;
Omega = 2*pi/365;

F1 = A(4)*(f - f_mean) + A(5)*(f - f_mean)^2 + A(6)*(f_mean - 150) + A(38)*(f_mean - 150)^2;
M = (A(7) + A(8)*P(1,3))*k;

if isT
    M = M + A(39)*exp(k);
else
    M = M + A(39)*k^2;
end

Q = A(2)*P(1,3) + A(3)*P(1,5) + A(37)*P(1,2);

AN1 = (A(9) + A(10)*P(1,3))*cos(Omega*(d - A(11)));
AN2 = (A(15)*P(1,2) + A(16)*P(1,4) + A(17)*P(1,6))*cos(Omega*(d - A(18)));

SAN1 = (A(12) + A(13)*P(1,3))*cos(2*Omega*(d - A(14)));
SAN2 = A(19)*P(1,2)*cos(2*Omega*(d - A(20)));

D = (A(21)*P(2,2) + A(22)*P(2,4) + A(23)*P(2,6) + ...
    (A(24)*P(2,2) + A(25)*P(2,3))*cos(Omega*(d - A(18))))*cos(omega*t) + ...
    (A(26)*P(2,2) + A(27)*P(2,4) + A(28)*P(2,6) + ...
    (A(29)*P(2,2) + A(30)*P(2,3))*cos(Omega*(d - A(18))))*sin(omega*t);

SD = (A(31)*P(3,3) + A(32)*P(3,4)*cos(Omega*(d - A(18))))*cos(2*omega*t) + ...
    (A(33)*P(3,3) + A(34)*P(3,4)*cos(Omega*(d - A(18))))*sin(2*omega*t);

TD = A(35)*P(4,4)*cos(3*omega*t) + A(36)*P(4,4)*sin(3*omega*t);

beta = 1;
if ~isH
    beta = beta + F1;
end

G = 1 + F1 + M + Q + beta*(AN1 + AN2 + SAN1 + SAN2 + D + SD + TD);

end