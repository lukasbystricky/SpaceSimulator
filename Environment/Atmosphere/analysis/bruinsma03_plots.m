%% Recreates plots in paper "The DTM-2000 empirical thermosphere model with
% new data assimilation and constraints at lower boundary: accuracy and 
% properties" by Bruinsma et al. 


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
% Temperature sensitivity to solar time and latitude
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
    
    data_lines = [min(min(T)), 750, 780, 800, 830, 860, 890, 920, 940, 970];
    contourf(t, lat*180/pi, T, data_lines);
    xlabel('solar time');
    ylabel('latitude  (^o)');
    title('200 km');
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
    
    data_lines = [min(min(T)), 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150];
    contourf(t, lat*180/pi, T, data_lines);
    xlabel('solar time');
    ylabel('latitude  (^o)');
    title('500 km');
    colorbar    
    
end

%% figure 9
% Density sensitivity to solar time and latitude
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
    
    data_lines = [min(min(rho)), [1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.7]*1e-11];
    contourf(t, lat*180/pi, rho, data_lines)
    
    xlabel('solar time');
    ylabel('latitude  (^o)');
    title('120 km');
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
    
    data_lines = [min(min(rho)), [2.2, 4.4, 6.6, 8.8, 11, 13, 15, 18, 20]*1e-17];
    contourf(t, lat*180/pi, rho, data_lines);
    
    xlabel('solar time');
    ylabel('latitude  (^o)');
    title('500 km');
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
    
    data_lines = [min(min(rho)), [1.4, 1.7, 1.9, 2.2, 2.4, 2.6, 2.9, 3.1, 3.4]*1e-18];
    contourf(t, lat*180/pi, rho, data_lines);
    
    xlabel('solar time');
    ylabel('latitude  (^o)');
    title('800 km');
    colorbar;

end
