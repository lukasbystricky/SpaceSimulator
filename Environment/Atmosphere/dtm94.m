function [t_inf, T, cH, cHe, cO, cN2] = dtm94(z, coeffs, P, f, f_mean, k, t, d)

T120 = 380;
T120p = 14.348;
R = 6.35677e3;

MU = 398600.4418;
g = MU/((R+120)^2);

BOLTZ = 1.38064852;

MASS_H = 1.6737236e2;
MASS_HE = 6.6464764e2;
MASS_O = 2.6567626e3;
MASS_N2 = 2*2.3258671e3;

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
a = (t_inf - T120)/t_inf;

T = t_inf - (t_inf - T120)*exp(-sigma*xi);

fH = (T120/T)^(1 + alphaH + gammaH)*exp(-sigma*gammaH*xi);
fHe = (T120/T)^(1 + alphaHe + gammaHe)*exp(-sigma*gammaHe*xi);
fO = (T120/T)^(1 + gammaO)*exp(-sigma*gammaO*xi);
fN2 = (T120/T)^(1 + gammaN2)*exp(-sigma*gammaN2*xi);

cH = cH_0*fH;
cHe = cHe_0*fHe;
cO = cO_0*fO;
cN2 = cN2_0*fN2;

% cH = cH_0*((1 - a)/(1 - a*exp(-sigma*xi)))^(1 + alphaH + gammaH)*exp(-sigma*xi*gammaH);
% cHe = cHe_0*((1 - a)/(1 - a*exp(-sigma*xi)))^(1 + alphaHe + gammaHe)*exp(-sigma*xi*gammaHe);
% cO = cO_0*((1 - a)/(1 - a*exp(-sigma*xi)))^(1 + gammaO)*exp(-sigma*xi*gammaO);
% cN2 = cN2_0*((1 - a)/(1 - a*exp(-sigma*xi)))^(1 + gammaN2)*exp(-sigma*xi*gammaN2);

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