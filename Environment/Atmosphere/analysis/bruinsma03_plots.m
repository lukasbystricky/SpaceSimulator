addpath ..
addpath ../coefficients
addpath ../../../Utilities/Math

close all

coeffs94 = dlmread('coeffs94.csv');

N = 30;

MASS_H = 1.6737236e-24;
MASS_HE = 6.6464764e-24;
MASS_O = 2.6567626e-23;
MASS_N2 = 2*2.3258671e-23;

FIGS_TO_PLOT = [6,9];

N = 50;
M = 100;
    
%% figure 6
if max(FIGS_TO_PLOT == 6)    
    
    lat = linspace(-pi/2, pi/2, N);
    
    d = 180;
    
    t = linspace(0,24,M);
    f = 150;
    k = 0;
    T = zeros(N,M);    
   
    z = 200;
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        
        for j = 1:M
            [ ~, T(i,j), ~, ~, ~, ~] = dtm94(...
                z, coeffs94, P, f, f, k, t(j), d);
        end
    end

    figure();
    subplot(1,2,1)
    contourf(t, lat*180/pi, T);
    colorbar
    
    z = 500;
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        
        for j = 1:M
            [ ~, T(i,j), ~, ~, ~, ~] = dtm94(...
                z, coeffs94, P, f, f, k, t(j), d);
        end
    end
    
    subplot(1,2,2)
    contourf(t, lat*180/pi, T);
    colorbar    
    
end

%% figure 9
if max(FIGS_TO_PLOT == 9)

    lat = linspace(-pi/2, pi/2, N);
    
    t = linspace(0,24,M);
    d = 15;
    
    f = 70;
    k = 0;   
    
    c_H = zeros(N,M);
    c_He = zeros(N,M);
    c_O = zeros(N,M);
    c_N2 = zeros(N,M);

    z = 120;
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        
        for j = 1:M
            [ ~, ~, c_H(i,j), c_He(i,j), c_O(i,j), c_N2(i,j)] = dtm94(...
                z, coeffs94, P, f, f, k, t(j), d);
        end
    end
    
    rho = (MASS_H*c_H + MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    figure();
    subplot(1,3,1);
    contourf(t, lat*180/pi, rho)
    colorbar;
    
    z = 500;
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        
        for j = 1:M
            [ ~, ~, c_H(i,j), c_He(i,j), c_O(i,j), c_N2(i,j)] = dtm94(...
                z, coeffs94, P, f, f, k, t(j), d);
        end
    end
    
    rho = (MASS_H*c_H + MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    subplot(1,3,2);
    contourf(t, lat*180/pi, rho)
    colorbar;
    
    z = 800;
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        
        for j = 1:M
            [ ~, ~, c_H(i,j), c_He(i,j), c_O(i,j), c_N2(i,j)] = dtm94(...
                z, coeffs94, P, f, f, k, t(j), d);
        end
    end
    
    rho = (MASS_H*c_H + MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    subplot(1,3,3);
    contourf(t, lat*180/pi, rho)
    colorbar;
    
end
