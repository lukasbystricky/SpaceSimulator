addpath ..
addpath ../coefficients
addpath ../../../Utilities/Math
addpath ../../../Utilities/Plots

close all

coeffs77 = dlmread('coeffs77.csv');

MASS_HE = 6.6464764e-24;
MASS_O = 2.6567626e-23;
MASS_N2 = 2*2.3258671e-23;
FIGS_TO_PLOT = [13];

%% figure 3
if max(FIGS_TO_PLOT == 3)
    
    N = 30;
    
    t = 8; %to be verified
    d = 264;
    k = 2;
    lat = 0;
    z = 400;   
    
    f = linspace(70,200,N);
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');    
    
    T = zeros(1,N);
    c_He = zeros(1,N);
    c_O = zeros(1,N);
    c_N2 = zeros(1,N);
    
    for i = 1:N
        [~, T(i), c_He(i), c_O(i), c_N2(i)] = dtm77(...
            z, coeffs77, P, f(i), f(i), k, t, d);
    end
    
    rho = (MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    figure();
    subplot(2,2,3);
    semilogy(f, rho, 'linewidth', 2);
    title('Density (g/cm^3)')
    xlabel('$f = \bar{f}$', 'Interpreter', 'LaTex');
    xlim([min(f), max(f)]);
    
    subplot(2,2,[2,4]);
    semilogy(f, 100*c_He, '-', 'linewidth', 2);
    hold
    semilogy(f, c_O, '-g', 'linewidth', 2);
    semilogy(f, c_N2, '-r', 'linewidth', 2);
    legend('100*He', 'O', 'N_2');
    title('Concentration (particles/cm^3)')
    xlabel('$f = \bar{f}$', 'Interpreter', 'LaTex');
    xlim([min(f), max(f)]);
    
    subplot(2,2,1);
    plot(f, T, 'linewidth', 2);
    title('Temperature (^oC)');
    xlabel('$f = \bar{f}$', 'Interpreter', 'LaTex');
    xlim([min(f), max(f)]);
    
end

%% figure 4
if max(FIGS_TO_PLOT == 4)
    
    N = 30;
    
    t = 8; %to be verified
    d = 264;
    k = 2;
    lat = 0;
    z = 400;
    f_mean = 130;
    
    f = linspace(100,180,N);
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    
    T = zeros(1,N);
    c_He = zeros(1,N);
    c_O = zeros(1,N);
    c_N2 = zeros(1,N);
    
    for i = 1:N
        [~, T(i), c_He(i), c_O(i), c_N2(i)] = dtm77(...
            z, coeffs77, P, f(i), f_mean, k, t, d);
    end
    
    rho = (MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    figure();
    subplot(2,2,3);
    semilogy(f - f_mean, rho, 'linewidth', 2);
    title('Density (g/cm^3)')
    xlabel('$f - \bar{f}$', 'Interpreter', 'LaTex');
    xlim([min(f - f_mean), max(f - f_mean)]);
    
    subplot(2,2,[2,4]);
    semilogy(f - f_mean, 100*c_He, '-', 'linewidth', 2);
    hold
    semilogy(f - f_mean, c_O, '-g', 'linewidth', 2);
    semilogy(f - f_mean, c_N2, '-r', 'linewidth', 2);
    legend('100*He', 'O', 'N_2');
    title('Concentration (particles/cm^3)')
    xlabel('$f - \bar{f}$', 'Interpreter', 'LaTex');
    xlim([min(f - f_mean), max(f - f_mean)]);
    
    subplot(2,2,1);
    plot(f - f_mean, T, 'linewidth', 2);
    title('Temperature (^oC)');
    xlabel('$f - \bar{f}$', 'Interpreter', 'LaTex');
    xlim([min(f - f_mean), max(f - f_mean)]);

end

%% figure 6
if max(FIGS_TO_PLOT == 6)

    N = 100;
    M = 100;
    lat = linspace(-pi/2, pi/2, N);
    
    t = 15;
    d =linspace(0, 365, M);
    
    % high solar activity
    f = 150;
    f_mean = 150;
    k = 2;
    
    z = 400;
    
    c_He = zeros(N,M);
    c_O = zeros(N,M);
    c_N2 = zeros(N,M);
    
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        
        for j = 1:M
            [ ~, ~, c_He(i,j), c_O(i,j), c_N2(i,j)] = dtm77(...
                z, coeffs77, P, f, f_mean, k, t, d(j));
        end
    end
    
    rho = (MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    figure();
    h = tight_subplot(4,2, [0.1, 0.01], [0.07, 0.05], 0.05);
    
    axes(h(2));
    contourf(d, lat*180/pi, rho);
    xlabel('Day');
    set(gca,'YTickLabel','')
    colorbar
    annotation('textbox', [0 0.9 1 0.1], 'String', 'Total density',...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    axes(h(4));
    contourf(d, lat*180/pi, c_N2);
    xlabel('Day');
    set(gca,'YTickLabel','')
    colorbar
    annotation('textbox', [0 0.665 1 0.1], 'String', 'N_2 concentration',...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    axes(h(6));
    contourf(d, lat*180/pi, c_O);
    xlabel('Day');
    set(gca,'YTickLabel','')
    colorbar
    annotation('textbox', [0 0.4 1 0.1], 'String', 'O concentration',...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    axes(h(8));
    contourf(d, lat*180/pi, c_He);
    xlabel('Day');
    set(gca,'YTickLabel','')
    colorbar
    annotation('textbox', [0 0.16 1 0.1], 'String', 'He concentration',...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    d = 180;
    t = linspace(0,24, M);
    
    c_He = zeros(N,M);
    c_O = zeros(N,M);
    c_N2 = zeros(N,M);
    
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        
        for j = 1:M
            [ ~, ~, c_He(i,j), c_O(i,j), c_N2(i,j)] = dtm77(...
                z, coeffs77, P, f, f_mean, k, t(j), d);
        end
    end
    
    rho = (MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    axes(h(1));
    contourf(t, lat*180/pi, rho);
    xlabel('Time');
    ylabel('Latitude');
    colorbar('westoutside');
    
    axes(h(3));
    contourf(t, lat*180/pi, c_N2);
    xlabel('Time');
    ylabel('Latitude');
    colorbar('westoutside');
    
    axes(h(5));
    contourf(t, lat*180/pi, c_O);
    xlabel('Time');
    ylabel('Latitude');
    colorbar('westoutside');
    
    axes(h(7));
    contourf(t, lat*180/pi, c_He);
    xlabel('Time');
    ylabel('Latitude');
    colorbar('westoutside');
    
end

%% figure 7
if max(FIGS_TO_PLOT == 7)
    
    N = 100;
    M = 100;
    lat = linspace(-pi/2, pi/2, N);
    
    t = 15;
    d =linspace(0, 365, M);
    
    % low solar activity
    f = 92;
    f_mean = 92;
    k = 1;
    
    z = 275;
    
    c_He = zeros(N,M);
    c_O = zeros(N,M);
    c_N2 = zeros(N,M);
    
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        
        for j = 1:M
            [ ~, ~, c_He(i,j), c_O(i,j), c_N2(i,j)] = dtm77(...
                z, coeffs77, P, f, f_mean, k, t, d(j));
        end
    end
    
    rho = (MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    figure();
    h = tight_subplot(4,2, [0.1, 0.01], [0.07, 0.05], 0.05);
    
    axes(h(2));
    contourf(d, lat*180/pi, rho);
    xlabel('Day');
    set(gca,'YTickLabel','')
    colorbar
    annotation('textbox', [0 0.9 1 0.1], 'String', 'Total density', 'EdgeColor',...
                'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    axes(h(4));
    contourf(d, lat*180/pi, c_N2);
    xlabel('Day');
    set(gca,'YTickLabel','')
    colorbar
    annotation('textbox', [0 0.665 1 0.1], 'String', 'N_2 concentration', 'EdgeColor',...
                'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    axes(h(6));
    contourf(d, lat*180/pi, c_O);
    xlabel('Day');
    set(gca,'YTickLabel','')
    colorbar
    annotation('textbox', [0 0.4 1 0.1], 'String', 'O concentration', 'EdgeColor',...
                'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    axes(h(8));
    contourf(d, lat*180/pi, c_He);
    xlabel('Day');
    set(gca,'YTickLabel','')
    colorbar
    annotation('textbox', [0 0.16 1 0.1], 'String', 'He concentration', 'EdgeColor',...
                'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    d = 180;
    t = linspace(0,24, M);
    
    c_He = zeros(N,M);
    c_O = zeros(N,M);
    c_N2 = zeros(N,M);
    
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        
        for j = 1:M
            [ ~, ~, c_He(i,j), c_O(i,j), c_N2(i,j)] = dtm77(...
                z, coeffs77, P, f, f_mean, k, t(j), d);
        end
    end
    
    rho = (MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    axes(h(1));
    contourf(t, lat*180/pi, rho);
    xlabel('Time');
    ylabel('Latitude');
    colorbar('westoutside');
    
    axes(h(3));
    contourf(t, lat*180/pi, c_N2);
    xlabel('Time');
    ylabel('Latitude');
    colorbar('westoutside');
    
    axes(h(5));
    contourf(t, lat*180/pi, c_O);
    xlabel('Time');
    ylabel('Latitude');
    colorbar('westoutside');
    
    axes(h(7));
    contourf(t, lat*180/pi, c_He);
    xlabel('Time');
    ylabel('Latitude');
    colorbar('westoutside');
    
end

%% figure 8
if max(FIGS_TO_PLOT == 8)
    
    d = 180;% to be verified
    k = 1;
    lat = 0;
    z = 400;
    f = 92;
    
    N = 30;
    t = linspace(0,24,N);
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    
    T = zeros(1,N);
    c_He = zeros(1,N);
    c_O = zeros(1,N);
    c_N2 = zeros(1,N);
    
    for i = 1:N
        [~, T(i), c_He(i), c_O(i), c_N2(i)] = dtm77(...
            z, coeffs77, P, f, f, k, t(i), d);
    end
    
    rho = (MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    rho_min = min(rho);
    rho_max = max(rho);
    c_He_min = min(c_He);
    c_He_max = max(c_He);
    c_O_min = min(c_O);
    c_O_max = max(c_O);
    c_N2_min = min(c_N2);
    c_N2_max = max(c_N2);
    
    figure();
    subplot(4,2,2);
    plot(t, c_N2);
    xlabel('Time');
    annotation('textbox', [0 0.9 1 0.1], 'String', 'N_2 concentration (particles/cm^3)', ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    subplot(4,2,4);
    plot(t, c_O);
    xlabel('Time');
    annotation('textbox', [0 0.675 1 0.1], 'String', 'O concentration (particles/cm^3)',...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    subplot(4,2,6);
    plot(t, c_He);
    xlabel('Time');
    annotation('textbox', [0 0.45 1 0.1], 'String', 'He concentration (particles/cm^3)',...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    subplot(4,2,8);
    plot(t, rho);
    xlabel('Time');
    annotation('textbox', [0 0.225 1 0.1], 'String', 'Total density (g/cm^3)', ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    % high solar activity    
    k = 2;
    f = 150;
    
    T = zeros(1,N);
    c_He = zeros(1,N);
    c_O = zeros(1,N);
    c_N2 = zeros(1,N);
    
    for i = 1:N
        [~, T(i), c_He(i), c_O(i), c_N2(i)] = dtm77(...
            z, coeffs77, P, f, f, k, t(i), d);
    end
    
    rho = (MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    rho_min = min(rho_min, min(rho));
    rho_max = max(rho_max, max(rho));
    c_He_min = min(c_He_min, min(c_He));
    c_He_max = max(c_He_max, max(c_He));
    c_O_min = min(c_O_min, min(c_O));
    c_O_max = max(c_O_max, max(c_O));
    c_N2_min = min(c_N2_min, min(c_N2));
    c_N2_max = max(c_N2_max, max(c_N2));
    
    subplot(4,2,1);
    plot(t, c_N2);
    xlabel('Time');
    ylim([c_N2_min, c_N2_max]);
    
    subplot(4,2,3);
    plot(t, c_O);
    xlabel('Time');
    ylim([c_O_min, c_O_max]);
    
    subplot(4,2,5);
    plot(t, c_He);
    xlabel('Time');
    ylim([c_He_min, c_He_max]);
    
    subplot(4,2,7);
    plot(t, rho);
    xlabel('Time');
    ylim([rho_min, rho_max]);
    
    subplot(4,2,2);
    ylim([c_N2_min, c_N2_max]);
    
    subplot(4,2,4);
    ylim([c_O_min, c_O_max]);
    
    subplot(4,2,6);
    ylim([c_He_min, c_He_max]);
    
    subplot(4,2,8)
    ylim([rho_min, rho_max]);
    
end

%% figures 9, 10, 11 and 12

if (max(FIGS_TO_PLOT == 9) || max(FIGS_TO_PLOT == 10) ||...
        max(FIGS_TO_PLOT == 11) || max(FIGS_TO_PLOT == 12))
    
    t = 15; %to be verified
    f = 150;
    k = 2;
    lat = 45*pi/180;
    z = 400;
    
    N = 30;
    d = linspace(0, 365, N);
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    
    T_inf_p45 = zeros(1,N);
    c_He_p45 = zeros(1,N);
    c_O_p45 = zeros(1,N);
    c_N2_p45 = zeros(1,N);
    
    for i = 1:N
        [T_inf_p45(i), ~, c_He_p45(i), c_O_p45(i), c_N2_p45(i)] = dtm77(...
            z, coeffs77, P, f, f, k, t, d(i));
    end
    
    rho_p45 = (MASS_HE*c_He_p45 + MASS_O*c_O_p45 + MASS_N2*c_N2_p45);
    
    lat = -45*pi/180;
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    
    T_inf_m45 = zeros(1,N);
    c_He_m45 = zeros(1,N);
    c_O_m45 = zeros(1,N);
    c_N2_m45 = zeros(1,N);
    
    for i = 1:N
        [T_inf_m45(i), ~, c_He_m45(i), c_O_m45(i), c_N2_m45(i)] = dtm77(...
            z, coeffs77, P, f, f, k, t, d(i));
    end
    
    rho_m45 = (MASS_HE*c_He_m45 + MASS_O*c_O_m45 + MASS_N2*c_N2_m45);
    
    figure();
    hold on
    plot(d, rho_p45);
    
    plot(d, rho_m45, 'r');
    
    legend('45^o latitude', '-45^o latitude');
    xlabel('day');
    
    title('Total density (g/cm^3)', 'Fontsize', 15);
    xlim([min(d), max(d)]);
    
    figure();
    subplot(1,3,3)
    hold on
    plot(d, c_He_p45);
    plot(d, c_He_m45);
    xlim([min(d), max(d)]);
    title('He concentration  (particles/cm^3)');
    xlabel('day');
    
    subplot(1,3,2)
    hold on
    plot(d, c_N2_p45);
    plot(d, c_N2_m45);
    xlim([min(d), max(d)]);
    title('N_2 concentration  (particles/cm^3)');
    xlabel('day');
    
    subplot(1,3,1)
    hold on
    plot(d, c_O_p45);
    plot(d, c_O_m45);
    xlim([min(d), max(d)]);
    title('O concentration  (particles/cm^3)');
    xlabel('day');
    
    figure();
    hold on
    plot(d, T_inf_p45);
    plot(d, T_inf_m45, 'r');
    
    legend('45^o latitude', '-45^o latitude');
    xlabel('day');
    
    title('Thermopause temperature (^oC)', 'Fontsize', 15);
    xlim([min(d), max(d)]);
    
end

%% figure 13
if max(FIGS_TO_PLOT == 13)
    
    t = 15; %to be verified
    f = 150;
    k = 2;
    d = 180;
    z = 1000;
    
    N = 30;
    lat = linspace(-pi/2, pi/2, N);
    
    c_He = zeros(1,N);
    
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        [~, ~, c_He(i), ~, ~] = dtm77(...
            z, coeffs77, P, f, f, k, t, d);
    end
    
    figure();
    
    subplot(1,2,1);
    semilogy(lat*180/pi, c_He);
    title('$f = \bar{f} = 150$, $k = 2$', 'interpreter', 'latex');
    hold on
    
    d = 360;
    c_He = zeros(1,N);
    
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        [~, ~, c_He(i), ~, ~] = dtm77(...
            z, coeffs77, P, f, f, k, t, d);
    end
    
    semilogy(lat*180/pi, c_He, 'r');
    
    legend('June 21', 'Dec 21', 'location', 'south');
    xlabel('latitude');
    ylabel('He concentration (particles/cm^3)');
    
    xlim([-90, 90]);
    ylim([2e4, 1.1e6]);
    
    d = 180;
    f = 92;
    k = 1;
    
    c_He = zeros(1,N);
    
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        [~, ~, c_He(i), ~, ~] = dtm77(...
            z, coeffs77, P, f, f, k, t, d);
    end
    
    subplot(1,2,2);
    
    semilogy(lat*180/pi, c_He);
    title('$f = \bar{f} = 92$, $k = 1$', 'interpreter', 'latex');
    hold on
    
    d = 360;
    c_He = zeros(1,N);
    
    for i = 1:N
        
        [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
        [~, ~, c_He(i), ~, ~] = dtm77(...
            z, coeffs77, P, f, f, k, t, d);
    end
    
    semilogy(lat*180/pi, c_He, 'r');
    
    xlabel('latitude');
    xlim([-90, 90]);
    ylim([2e4, 1.1e6]);

end

