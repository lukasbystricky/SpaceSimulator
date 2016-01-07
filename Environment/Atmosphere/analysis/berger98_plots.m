%% Recreates plots in paper "Improvement of the empirical thermospheric 
% model DTM: DTM94 - a comparative review of various temporal variations
% and prospects in space geodesy applications" by Berger et al. 

addpath ..
addpath ../coefficients
addpath ../../../Utilities/Math

close all

coeffs77 = dlmread('coeffs77.csv');
coeffs94 = dlmread('coeffs94.csv');

N = 30;

MASS_H = 1.6737236e-24;
MASS_HE = 6.6464764e-24;
MASS_O = 2.6567626e-23;
MASS_N2 = 2*2.3258671e-23;

FIGS_TO_PLOT = [9];

%% figure 2
% Geomagnetic sensitity for density at various altitudes
if max(FIGS_TO_PLOT == 2)
    
    d = 80;
    f = 150;
    lat = 45*pi/180;
    t = 9;
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
   
    k = linspace(0, 8.5, N);
    
    rho_800_77 = zeros(1,N);
    rho_800_94 = zeros(1,N);
    rho_1400_77 = zeros(1,N);
    rho_1400_94 = zeros(1,N);
    rho_300_77 = zeros(1,N);
    rho_300_94 = zeros(1,N);
    rho_500_77 = zeros(1,N);
    rho_500_94 = zeros(1,N);
    
    z = 800;
    for i = 1:N
        [~, ~, c_He_77, c_O_77, c_N2_77] = dtm77(...
            z, coeffs77, P, f, f, k(i), t, d);
        [~, ~, c_H_94, c_He_94, c_O_94, c_N2_94] = dtm94(...
            z, coeffs94, P, f, f, k(i), t, d);
        
        rho_800_77(i) = MASS_HE*c_He_77 + MASS_O*c_O_77 + MASS_N2*c_N2_77;
        rho_800_94(i) = MASS_H*c_H_94 + MASS_HE*c_He_94 + MASS_O*c_O_94 + MASS_N2*c_N2_94;
        
    end
    
    z = 1400;
    for i = 1:N
         [~, ~, c_He_77, c_O_77, c_N2_77] = dtm77(...
            z, coeffs77, P, f, f, k(i), t, d);
        [~, ~, c_H_94, c_He_94, c_O_94, c_N2_94] = dtm94(...
            z, coeffs94, P, f, f, k(i), t, d);
        
        rho_1400_77(i) = MASS_HE*c_He_77 + MASS_O*c_O_77 + MASS_N2*c_N2_77;
        rho_1400_94(i) = MASS_H*c_H_94 + MASS_HE*c_He_94 + MASS_O*c_O_94 + MASS_N2*c_N2_94;
    end
    
    z = 300;
    for i = 1:N
         [~, ~, c_He_77, c_O_77, c_N2_77] = dtm77(...
            z, coeffs77, P, f, f, k(i), t, d);
        [~, ~, c_H_94, c_He_94, c_O_94, c_N2_94] = dtm94(...
            z, coeffs94, P, f, f, k(i), t, d);
        
        rho_300_77(i) = MASS_HE*c_He_77 + MASS_O*c_O_77 + MASS_N2*c_N2_77;
        rho_300_94(i) = MASS_H*c_H_94 + MASS_HE*c_He_94 + MASS_O*c_O_94 + MASS_N2*c_N2_94;
    end
    
    z = 500;
    for i = 1:N
         [~, ~, c_He_77, c_O_77, c_N2_77] = dtm77(...
            z, coeffs77, P, f, f, k(i), t, d);
        [~, ~, c_H_94, c_He_94, c_O_94, c_N2_94] = dtm94(...
            z, coeffs94, P, f, f, k(i), t, d);
        
        rho_500_77(i) = MASS_HE*c_He_77 + MASS_O*c_O_77 + MASS_N2*c_N2_77;
        rho_500_94(i) = MASS_H*c_H_94 + MASS_HE*c_He_94 + MASS_O*c_O_94 + MASS_N2*c_N2_94;
    end
    
    figure();
    
    subplot(2,2,1)
    plot(k, rho_800_77);
    hold on
    plot(k, rho_800_94, 'r');
    legend('DTM 77', 'DTM 94','Location','northwest');
    xlabel('Geomagnetic index');
    ylabel('g/cm^3');
    title('800 km');
    
    subplot(2,2,2)
    plot(k, rho_1400_77);
    hold on
    plot(k, rho_1400_94, 'r');
    xlabel('Geomagnetic index');
    ylabel('g/cm^3');
    title('1400 km');
    
    subplot(2,2,3)
    plot(k, rho_300_77);
    hold on
    plot(k, rho_300_94, 'r');
    xlabel('Geomagnetic index');
    ylabel('g/cm^3');
    title('300 km');
    
    subplot(2,2,4)
    plot(k, rho_500_77);
    hold on
    plot(k, rho_500_94, 'r');
    xlabel('Geomagnetic index');
    ylabel('g/cm^3');
    title('500 km');   
    
end

%% figure 3
% Sensitivity of concentrations to geomagnetic index
if max(FIGS_TO_PLOT == 3)
    d = 80;
    f = 150;
    t = 9;    
    lat = 45*pi/180;
    
    k = linspace(0,8.5,N);
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    
    c_He_77 = zeros(1,N);
    c_O_77 = zeros(1,N);
    c_N2_77 = zeros(1,N);
    c_He_94 = zeros(1,N);
    c_O_94 = zeros(1,N);
    c_N2_94 = zeros(1,N);
    c_He_77_80 = zeros(1,N);
    c_He_94_80 = zeros(1,N);
    
    z = 500;
    for i = 1:N
        [~, ~, ~, c_O_77(i), ~] = dtm77(...
            z, coeffs77, P, f, f, k(i), t, d);
        [~, ~, ~, ~, c_O_94(i), ~] = dtm94(...
            z, coeffs94, P, f, f, k(i), t, d);
    end
    
    z = 300;
    for i = 1:N
        [~, ~, ~, ~, c_N2_77(i)] = dtm77(...
            z, coeffs77, P, f, f, k(i), t, d);
        [~, ~, ~, ~, ~, c_N2_94(i)] = dtm94(...
            z, coeffs94, P, f, f, k(i), t, d);
    end
    
    lat = 0;
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    z = 1400;
    for i = 1:N
        [~, ~, c_He_77(i), ~, ~] = dtm77(...
            z, coeffs77, P, f, f, k(i), t, d);
        [~, ~, ~, c_He_94(i), ~, ~] = dtm94(...
            z, coeffs94, P, f, f, k(i), t, d);
    end   
    
    lat = 80*pi/180;
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    for i = 1:N
        [~, ~, c_He_77_80(i), ~, ~] = dtm77(...
            z, coeffs77, P, f, f, k(i), t, d);
        [~, ~, ~, c_He_94_80(i), ~, ~] = dtm94(...
            z, coeffs94, P, f, f, k(i), t, d);
    end
    
    figure();
    
    subplot(2,2,1);
    plot(k, c_O_77);
    hold on
    plot(k, c_O_94, 'r');
    legend('DTM 77', 'DTM 94','Location','northwest');
    title('O concentration (500 km)');
    xlabel('Geomagnetic index');
    ylabel('Particles per cm^3');
    
    subplot(2,2,2);
    plot(k, c_N2_77);
    hold on
    plot(k, c_N2_94, 'r');
    title('N_2 concentration (300 km)');
    xlabel('Geomagnetic index');
    ylabel('Particles per cm^3');
    
    subplot(2,2,3);
    plot(k, c_He_77)
    hold on
    plot(k, c_He_94, 'r');
    title('He concentration (1400 km, equator)');
    xlabel('Geomagnetic index');
    ylabel('Particles per cm^3');
    
    subplot(2,2,4);
    plot(k, c_He_77_80)
    hold on
    plot(k, c_He_94_80, 'r');
    title('He concentration (1400 km, 80^o latitude)');
    xlabel('Geomagnetic index');
    ylabel('Particles per cm^3');
    
end

%% figure 4
% Sensitivity of the thermopause temperature to geomagnetic index and
% solar flux
if max(FIGS_TO_PLOT == 4)

    d = 80;
    f = 150;
    lat = 45*pi/180;
    t = 9;
    z = 500;
    k = linspace(0, 8.5,N);
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    
    T_inf_77_geo = zeros(1,N);
    T_inf_94_geo = zeros(1,N);
    T_inf_77_flux = zeros(1,N);
    T_inf_94_flux = zeros(1,N);
    
   
    for i = 1:N
        [T_inf_77_geo(i), ~, ~, ~, ~] = dtm77(...
            z, coeffs77, P, f, f, k(i), t, d);
        [T_inf_94_geo(i), ~, ~, ~, ~, ~] = dtm94(...
            z, coeffs94, P, f, f, k(i), t, d);
    end
    
    figure();
    subplot(1,2,1)
    plot(k, T_inf_77_geo);
    hold on
    plot(k, T_inf_94_geo, 'r');
    
    legend('DTM 77', 'DTM 94','Location','northwest');
    xlabel('Geomagnetic index');

    
    k = 3;
    f = linspace(50, 250, N);
    
    for i = 1:N
        [T_inf_77_flux(i), ~, ~, ~, ~] = dtm77(...
            z, coeffs77, P, f(i), f(i), k, t, d);
        [T_inf_94_flux(i), ~, ~, ~, ~, ~] = dtm94(...
            z, coeffs94, P, f(i), f(i), k, t, d);
    end
    
    subplot(1,2,2)
    plot(f, T_inf_77_flux);
    hold on
    plot(f, T_inf_94_flux, 'r');
    
    xlabel('Mean solar flux $f = \bar{f}$', 'interpreter' ,'latex');   
    
    annotation('textbox', [0 0.9 1 0.1], 'String', 'Thermopause temperature (^oK)',...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
end

%% figure 5
% Density sensitivity to f = f_mean
if max(FIGS_TO_PLOT == 5)
    
    d = 80;
    k = 3;
    lat = 45*pi/180;
    t = 9;
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
   
    f = linspace(50, 250, N);
    
    rho_800_77 = zeros(1,N);
    rho_800_94 = zeros(1,N);
    rho_1400_77 = zeros(1,N);
    rho_1400_94 = zeros(1,N);
    rho_300_77 = zeros(1,N);
    rho_300_94 = zeros(1,N);
    rho_500_77 = zeros(1,N);
    rho_500_94 = zeros(1,N);
    
    z = 800;
    for i = 1:N
        [~, ~, c_He_77, c_O_77, c_N2_77] = dtm77(...
            z, coeffs77, P, f(i), f(i), k, t, d);
        [~, ~, c_H_94, c_He_94, c_O_94, c_N2_94] = dtm94(...
            z, coeffs94, P, f(i), f(i), k, t, d);
        
        rho_800_77(i) = MASS_HE*c_He_77 + MASS_O*c_O_77 + MASS_N2*c_N2_77;
        rho_800_94(i) = MASS_H*c_H_94 + MASS_HE*c_He_94 + MASS_O*c_O_94 + MASS_N2*c_N2_94;
        
    end
    
    z = 1400;
    for i = 1:N
        [~, ~, c_He_77, c_O_77, c_N2_77] = dtm77(...
            z, coeffs77, P, f(i), f(i), k, t, d);
        [~, ~, c_H_94, c_He_94, c_O_94, c_N2_94] = dtm94(...
            z, coeffs94, P, f(i), f(i), k, t, d);
        
        rho_1400_77(i) = MASS_HE*c_He_77 + MASS_O*c_O_77 + MASS_N2*c_N2_77;
        rho_1400_94(i) = MASS_H*c_H_94 + MASS_HE*c_He_94 + MASS_O*c_O_94 + MASS_N2*c_N2_94;
    end
    
    z = 300;
    for i = 1:N
        [~, ~, c_He_77, c_O_77, c_N2_77] = dtm77(...
            z, coeffs77, P, f(i), f(i), k, t, d);
        [~, ~, c_H_94, c_He_94, c_O_94, c_N2_94] = dtm94(...
            z, coeffs94, P, f(i), f(i), k, t, d);
        
        rho_300_77(i) = MASS_HE*c_He_77 + MASS_O*c_O_77 + MASS_N2*c_N2_77;
        rho_300_94(i) = MASS_H*c_H_94 + MASS_HE*c_He_94 + MASS_O*c_O_94 + MASS_N2*c_N2_94;
    end
    
    z = 500;
    for i = 1:N
        [~, ~, c_He_77, c_O_77, c_N2_77] = dtm77(...
            z, coeffs77, P, f(i), f(i), k, t, d);
        [~, ~, c_H_94, c_He_94, c_O_94, c_N2_94] = dtm94(...
            z, coeffs94, P, f(i), f(i), k, t, d);
        
        rho_500_77(i) = MASS_HE*c_He_77 + MASS_O*c_O_77 + MASS_N2*c_N2_77;
        rho_500_94(i) = MASS_H*c_H_94 + MASS_HE*c_He_94 + MASS_O*c_O_94 + MASS_N2*c_N2_94;
    end
    
    figure();
    
    subplot(2,2,1)
    plot(f, rho_800_77);
    hold on
    plot(f, rho_800_94, 'r');
    legend('DTM 77', 'DTM 94','Location','northwest');
    xlabel('Mean solar flux');
    ylabel('g/cm^3');
    title('800 km');
    
    subplot(2,2,2)
    plot(f, rho_1400_77);
    hold on
    plot(f, rho_1400_94, 'r');
    xlabel('Mean solar flux');
    ylabel('g/cm^3');
    title('1400 km');
    
    subplot(2,2,3)
    plot(f, rho_300_77);
    hold on
    plot(f, rho_300_94, 'r');
    xlabel('Mean solar flux');
    ylabel('g/cm^3');
    title('300 km');
    
    subplot(2,2,4)
    plot(f, rho_500_77);
    hold on
    plot(f, rho_500_94, 'r');
    xlabel('Mean solar flux');
    ylabel('g/cm^3');
    title('500 km');    
end

%% figure 6
% Density and thermopause temperature sensitivity to f - f_mean
if max(FIGS_TO_PLOT == 6)
    
    d = 80;
    k = 3;
    lat = 45*pi/180;
    t = 9;
    f_mean = 100;
    z = 500;
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
   
    f_day = linspace(75, 175, N);
    
    rho_77 = zeros(1,N);
    rho_94 = zeros(1,N);
    T_inf_77 = zeros(1,N);
    T_inf_94 = zeros(1,N);

    for i = 1:N
        [T_inf_77(i), ~, c_He_77, c_O_77, c_N2_77] = dtm77(...
            z, coeffs77, P, f_day(i), f_mean, k, t, d);
        [T_inf_94(i), ~, c_H_94, c_He_94, c_O_94, c_N2_94] = dtm94(...
            z, coeffs94, P, f_day(i), f_mean, k, t, d);
        
        rho_77(i) = MASS_HE*c_He_77 + MASS_O*c_O_77 + MASS_N2*c_N2_77;
        rho_94(i) = MASS_H*c_H_94 + MASS_HE*c_He_94 + MASS_O*c_O_94 + MASS_N2*c_N2_94;
        
    end
    
    figure();
    
    subplot(2,2,1)
    plot(f_day - f_mean, rho_77);
    hold on
    plot(f_day - f_mean, rho_94, 'r');
    legend('DTM 77', 'DTM 94','Location','northwest');
    xlabel('$f - \bar{f}$', 'Interpreter', 'LaTex');
    ylabel('g/cm^3');
    title('Total density, $\bar{f} = 100$', 'Interpreter', 'LaTex');
    
    subplot(2,2,3)
    plot(f_day - f_mean, T_inf_77);
    hold on
    plot(f_day - f_mean, T_inf_94, 'r');
    xlabel('$f - \bar{f}$', 'Interpreter', 'LaTex');
    ylabel('^o K');
    title('Thermopause temperature, $\bar{f} = 100$', 'Interpreter', 'LaTex');
    
    f_mean = 200;
    f_day = linspace(125, 275, N);
    
    rho_77 = zeros(1,N);
    rho_94 = zeros(1,N);
    T_inf_77 = zeros(1,N);
    T_inf_94 = zeros(1,N);
    
    for i = 1:N
        [T_inf_77(i), ~, c_He_77, c_O_77, c_N2_77] = dtm77(...
            z, coeffs77, P, f_day(i), f_mean, k, t, d);
        [T_inf_94(i), ~, c_H_94, c_He_94, c_O_94, c_N2_94] = dtm94(...
            z, coeffs94, P, f_day(i), f_mean, k, t, d);
        
        rho_77(i) = MASS_HE*c_He_77 + MASS_O*c_O_77 + MASS_N2*c_N2_77;
        rho_94(i) = MASS_H*c_H_94 + MASS_HE*c_He_94 + MASS_O*c_O_94 + MASS_N2*c_N2_94;
        
    end
    
    subplot(2,2,2)
    plot(f_day - f_mean, rho_77);
    hold on
    plot(f_day - f_mean, rho_94, 'r');
    xlabel('$f - \bar{f}$',  'Interpreter', 'LaTex');
    ylabel('g/cm^3');
    title('Total density, $\bar{f} = 200$',  'Interpreter', 'LaTex');
    
    subplot(2,2,4)
    plot(f_day - f_mean, T_inf_77);
    hold on
    plot(f_day - f_mean, T_inf_94, 'r');
    xlabel('$f - \bar{f}$', 'interpreter', 'latex');
    ylabel('^o K');
    title('Thermopause temperature, $\bar{f} = 200$',  'Interpreter', 'LaTex');
    
end

%% figure 7
% Solar time sensitivity for concentrations and thermopause temperature
if max(FIGS_TO_PLOT == 7)

    d = 80;
    f = 150;
    k = 3;
    lat = 45*pi/180;
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    
    t = linspace(0, 24,N);
    
    T_inf_77 = zeros(1,N);
    c_He_77 = zeros(1,N);
    c_O_77 = zeros(1,N);
    c_N2_77 = zeros(1,N);
    T_inf_94 = zeros(1,N);
    c_He_94 = zeros(1,N);
    c_O_94 = zeros(1,N);
    c_N2_94 = zeros(1,N);
    
    z = 500;
    for i = 1:N
        [T_inf_77(i), ~, ~, c_O_77(i), ~] = dtm77(...
            z, coeffs77, P, f, f, k, t(i), d);
        [T_inf_94(i), ~, ~, ~, c_O_94(i), ~] = dtm94(...
            z, coeffs94, P, f, f, k, t(i), d);
    end
    
    z = 1400;
    for i = 1:N
        [~, ~, c_He_77(i), ~, ~] = dtm77(...
            z, coeffs77, P, f, f, k, t(i), d);
        [~, ~, ~, c_He_94(i), ~, ~] = dtm94(...
            z, coeffs94, P, f, f, k, t(i), d);
    end
    
    z = 300;
    for i = 1:N
        [~, ~, ~, ~, c_N2_77(i)] = dtm77(...
            z, coeffs77, P, f, f, k, t(i), d);
        [~, ~, ~, ~, ~, c_N2_94(i)] = dtm94(...
            z, coeffs94, P, f, f, k, t(i), d);
    end
    
    figure();
    
    subplot(2,2,1);
    plot(t, c_O_77);
    hold on
    plot(t, c_O_94, 'r');
    xlim([0,24]);
    legend('DTM 77', 'DTM 94');
    title('O concentration (500 km)');
    xlabel('solar time');
    ylabel('particles per cm^3');
    
    subplot(2,2,2);
    plot(t, c_He_77);
    hold on
    plot(t, c_He_94, 'r');
    xlim([0,24]);
    ylim([0, 3e5]);
    title('He concentration (1400 km)');
    xlabel('solar time');
    ylabel('particles per cm^3');
    
    subplot(2,2,3);
    plot(t, T_inf_77)
    hold on
    plot(t, T_inf_94, 'r');
    xlim([0,24]);
    title('Thermopause temperature');
    xlabel('solar time');
    ylabel('degree Kelvin');
    
    subplot(2,2,4);
    plot(t, c_N2_77)
    hold on
    plot(t, c_N2_94, 'r');
    xlim([0,24]);
    title('N_2 concentration (300 km)');
    xlabel('solar time');
    ylabel('particles per cm^3');
    
end

%% figure 8
% Day of year sensitiviy for concentrations and theropause temperature
if max(FIGS_TO_PLOT == 8)

    t = 9;
    f = 150;
    k = 3;
    lat = 45*pi/180;
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    
    d = linspace(0, 365 ,N);
    
    T_inf_77 = zeros(1,N);
    c_He_77 = zeros(1,N);
    c_O_77 = zeros(1,N);
    c_N2_77 = zeros(1,N);
    T_inf_94 = zeros(1,N);
    c_He_94 = zeros(1,N);
    c_O_94 = zeros(1,N);
    c_N2_94 = zeros(1,N);
    
    z = 500;
    for i = 1:N
        [T_inf_77(i), ~, ~, c_O_77(i), ~] = dtm77(...
            z, coeffs77, P, f, f, k, t, d(i));
        [T_inf_94(i), ~, ~, ~, c_O_94(i), ~] = dtm94(...
            z, coeffs94, P, f, f, k, t, d(i));
    end
    
    z = 1400;
    for i = 1:N
        [~, ~, c_He_77(i), ~, ~] = dtm77(...
            z, coeffs77, P, f, f, k, t, d(i));
        [~, ~, ~, c_He_94(i), ~, ~] = dtm94(...
            z, coeffs94, P, f, f, k, t, d(i));
    end
    
    z = 300;
    for i = 1:N
        [~, ~, ~, ~, c_N2_77(i)] = dtm77(...
            z, coeffs77, P, f, f, k, t, d(i));
        [~, ~, ~, ~, ~, c_N2_94(i)] = dtm94(...
            z, coeffs94, P, f, f, k, t, d(i));
    end
    
    figure();
    
    subplot(2,2,1);
    plot(d, c_O_77);
    hold on
    plot(d, c_O_94, 'r');
    xlim([0,365]);
    legend('DTM 77', 'DTM 94');
    title('O concentration (500 km)');
    xlabel('day of year');
    ylabel('particles per cm^3');
    
    subplot(2,2,2);
    plot(d, c_He_77);
    hold on
    plot(d, c_He_94, 'r');
    xlim([0,365]);
    title('He concentration (1400 km)');
    xlabel('day of year');
    ylabel('particles per cm^3');
    
    subplot(2,2,3);
    plot(d, T_inf_77)
    hold on
    plot(d, T_inf_94, 'r');
    xlim([0,365]);
    title('Thermopause temperature');
    xlabel('day of year');
    ylabel('degree Kelvin');
    
    subplot(2,2,4);
    plot(d, c_N2_77)
    hold on
    plot(d, c_N2_94, 'r');
    xlim([0,365]);
    title('N_2 concentration (300 km)');
    xlabel('day of year');
    ylabel('particles per cm^3');
    
end

%% figure 9
% Latitude sensitivity for concentrations and thermopause temperature
if max(FIGS_TO_PLOT == 9)

    d = 80;
    t = linspace(0, 24, N);
    f = 150;
    k = 3;
    lat = linspace(-pi/2, pi/2, N);
    
 
    
    T_inf_77 = zeros(1,N);
    c_He_77 = zeros(1,N);
    c_O_77 = zeros(1,N);
    c_N2_77 = zeros(1,N);
    T_inf_94 = zeros(1,N);
    c_He_94 = zeros(1,N);
    c_O_94 = zeros(1,N);
    c_N2_94 = zeros(1,N);
    
    z = 500;
    for i = 1:N   
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        
        for j = 1:N
               [t_part_77, ~, ~, c_part_77, ~] = dtm77(z, coeffs77, P, f, f, k, t(j), d);
               [t_part_94, ~, ~, ~, c_part_94, ~] = dtm94(z, coeffs94, P, f, f, k, t(j), d);
               
               T_inf_77(i) = T_inf_77(i) + t_part_77;
               c_O_77(i) = c_O_77(i) + c_part_77;
               
               T_inf_94(i) = T_inf_94(i) + t_part_94;
               c_O_94(i) = c_O_94(i) + c_part_94;
        end
    end
    
    z = 1400;
    for i = 1:N
       [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        
        for j = 1:N
               [~, ~, c_part_77, ~,  ~] = dtm77(z, coeffs77, P, f, f, k, t(j), d);
               [~, ~, ~, c_part_94,~,  ~] = dtm94(z, coeffs94, P, f, f, k, t(j), d);
               
               c_He_77(i) = c_He_77(i) + c_part_77;
               c_He_94(i) = c_He_94(i) + c_part_94;
        end
    end
    
    z = 300;
    for i = 1:N
       [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        
        for j = 1:N
               [~, ~, ~,  ~, c_part_77] = dtm77(z, coeffs77, P, f, f, k, t(j), d);
               [~, ~, ~,~,  ~, c_part_94] = dtm94(z, coeffs94, P, f, f, k, t(j), d);
               
               c_N2_77(i) = c_N2_77(i) + c_part_77;
               c_N2_94(i) = c_N2_94(i) + c_part_94;
        end
    end
    
    figure();
    
    subplot(2,2,1);
    plot(lat*180/pi, c_O_77/N);
    hold on
    plot(lat*180/pi, c_O_94/N, 'r');
    legend('DTM 77', 'DTM 94','Location','northwest');
    xlim([-90,90]);
    ylim([0, 7e7]);
    title('O concentration (500 km)');
    xlabel('latitude');
    ylabel('particles per cm^3');
    
    subplot(2,2,2);
    plot(lat*180/pi, c_He_77/N);
    hold on
    plot(lat*180/pi, c_He_94/N, 'r');
    xlim([-90,90]);
    ylim([0, 3e5]);
    title('He concentration (1400 km)');
    xlabel('latitude');
    ylabel('particles per cm^3');
    
    subplot(2,2,3);
    plot(lat*180/pi, T_inf_77/N)
    hold on
    plot(lat*180/pi, T_inf_94/N, 'r');
    xlim([-90,90]);
    ylim([900, 1350]);
    title('Thermopause temperature');
    xlabel('latitude');
    ylabel('degree Kelvin');
    
    subplot(2,2,4);
    plot(lat*180/pi, c_N2_77/N)
    hold on
    plot(lat*180/pi, c_N2_94/N, 'r');
    xlim([-90,90]);
    ylim([0, 3e8]);
    title('N_2 concentration (300 km)');
    xlabel('latitude');
    ylabel('particles per cm^3');
end