addpath ..
addpath ../coefficients
addpath ../../../Utilities/Math
addpath ../../../Utilities/Plots

close all

coeffs77 = dlmread('coeffs77_update.csv');

MASS_HE = 6.6464764e-24;
MASS_O = 2.6567626e-23;
MASS_N2 = 2*2.3258671e-23;
FIGS_TO_PLOT = [5];

%% figure 3
% Sensitivity of density, temperature and concentrations to f = f_mean
if max(FIGS_TO_PLOT == 3)
    
    N = 30;
    
    t = linspace(0,24,N);
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
        for j = 1:N
            [~, T_tmp, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
                z, coeffs77, P, f(i), f(i), k, t(j), d);
            
            T(i) = T(i) + T_tmp/N;
            c_He(i) = c_He(i) + c_He_tmp/N;
            c_O(i) = c_O(i) + c_O_tmp/N;
            c_N2(i) = c_N2(i) + c_N2_tmp/N;
        end
    end
    
    rho = (MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    figure();
    subplot(2,2,3);
    semilogy(f, rho, 'linewidth', 2);
    title('Density (g/cm^3)')
    xlabel('$f = \bar{f}$ ($10^{-22}$W m$^{-2}$Hz$^{-1}$)', 'Interpreter', 'LaTex');
    xlim([min(f), max(f)]);
    
    subplot(2,2,[2,4]);
    semilogy(f, 100*c_He, '-', 'linewidth', 2);
    hold
    semilogy(f, c_O, '-g', 'linewidth', 2);
    semilogy(f, c_N2, '-r', 'linewidth', 2);
    legend('100*He', 'O', 'N_2');
    title('Concentration (particles/cm^3)')
    xlabel('$f = \bar{f}$ ($10^{-22}$W m$^{-2}$Hz$^{-1}$)', 'Interpreter', 'LaTex');
    xlim([min(f), max(f)]);
    
    subplot(2,2,1);
    plot(f, T, 'linewidth', 2);
    title('Temperature (^oK)');
    xlabel('$f = \bar{f}$ ($10^{-22}$W m$^{-2}$Hz$^{-1}$)', 'Interpreter', 'LaTex');
    xlim([min(f), max(f)]);
    
end


%% figure 4
% Sensitivity of density, temperature and concentrations to f - f_mean
if max(FIGS_TO_PLOT == 4)
    
    N = 30;
    
    t = linspace(0,24,N);
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
        for j = 1:N
            [~, T_tmp, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
                    z, coeffs77, P, f(i), f_mean, k, t(j), d);
                
            T(i) = T(i) + T_tmp/N;
            c_He(i) = c_He(i) + c_He_tmp/N;
            c_O(i) = c_O(i) + c_O_tmp/N;
            c_N2(i) = c_N2(i) + c_N2_tmp/N;    
        end
    end
    
    rho = (MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    figure();
    subplot(2,2,3);
    semilogy(f - f_mean, rho, 'linewidth', 2);
    title('Density (g/cm^3)')
    xlabel('$f - \bar{f}$($10^{-22}$W m$^{-2}$Hz$^{-1}$)', 'Interpreter', 'LaTex');
    xlim([min(f - f_mean), max(f - f_mean)]);
    
    subplot(2,2,[2,4]);
    semilogy(f - f_mean, 100*c_He, '-', 'linewidth', 2);
    hold
    semilogy(f - f_mean, c_O, '-g', 'linewidth', 2);
    semilogy(f - f_mean, c_N2, '-r', 'linewidth', 2);
    legend('100*He', 'O', 'N_2');
    title('Concentration (particles/cm^3)')
    xlabel('$f - \bar{f}$ ($10^{-22}$W m$^{-2}$Hz$^{-1}$)', 'Interpreter', 'LaTex');
    xlim([min(f - f_mean), max(f - f_mean)]);
    
    subplot(2,2,1);
    plot(f - f_mean, T, 'linewidth', 2);
    title('Temperature (^oK)');
    xlabel('$f - \bar{f}$ ($10^{-22}$W m$^{-2}$Hz$^{-1}$)', 'Interpreter', 'LaTex');
    xlim([min(f - f_mean), max(f - f_mean)]);

end

%% figure 5
if max(FIGS_TO_PLOT == 5)
   N = 30;
   lat = linspace(0, pi/2, N);
   
   d = 80;
   f_1 = 114;
   f_2 = 112;
   
   k_1 = 0;
   k_2 = 5;
   
   t = linspace(0,24,N);
   
   z = 200;
  
   c_He_1 = zeros(1,N);
   c_O_1 = zeros(1,N);
   c_N2_1 = zeros(1,N);

   c_He_2 = zeros(1,N);
   c_O_2 = zeros(1,N);
   c_N2_2 = zeros(1,N);
   
   for i = 1:N
       [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
       
       for j = 1:N
           [~, ~, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
               z, coeffs77, P, f_1, f_1, k_1, t(j), d);
           
           c_He_1(i) = c_He_1(i) + c_He_tmp/N;
           c_O_1(i) = c_O_1(i) + c_O_tmp/N;
           c_N2_1(i) = c_N2_1(i) + c_N2_tmp/N;
           
           [~, ~, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
               z, coeffs77, P, f_2, f_2, k_2, t(j), d);
         
          
           c_He_2(i) = c_He_2(i) + c_He_tmp/N;
           c_O_2(i) = c_O_2(i) + c_O_tmp/N;
           c_N2_2(i) = c_N2_2(i) + c_N2_tmp/N;
                 
       end
   end
   
   rho_1 = (MASS_HE*c_He_1 + MASS_O*c_O_1 + MASS_N2*c_N2_1);
   rho_2 = (MASS_HE*c_He_2 + MASS_O*c_O_2 + MASS_N2*c_N2_2); 
   
   figure();
   subplot(2,2,1)
   semilogy(180*lat/pi, c_N2_2./c_N2_1);
   
   subplot(2,2,2)
   semilogy(180*lat/pi, c_O_2./c_O_1);
   
   subplot(2,2,3)
   semilogy(180*lat/pi, rho_2./rho_1);
   
   subplot(2,2,4)
   semilogy(180*lat/pi, c_He_2./c_He_1);

   z = 400;
   
   c_He_1 = zeros(1,N);
   c_O_1 = zeros(1,N);
   c_N2_1 = zeros(1,N);

   c_He_2 = zeros(1,N);
   c_O_2 = zeros(1,N);
   c_N2_2 = zeros(1,N);
   
   for i = 1:N
       [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
       
       for j = 1:N
           [~, ~, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
               z, coeffs77, P, f_1, f_1, k_1, t(j), d);
           
           c_He_1(i) = c_He_1(i) + c_He_tmp/N;
           c_O_1(i) = c_O_1(i) + c_O_tmp/N;
           c_N2_1(i) = c_N2_1(i) + c_N2_tmp/N;
           
           [~, ~, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
               z, coeffs77, P, f_2, f_2, k_2, t(j), d);
         
          
           c_He_2(i) = c_He_2(i) + c_He_tmp/N;
           c_O_2(i) = c_O_2(i) + c_O_tmp/N;
           c_N2_2(i) = c_N2_2(i) + c_N2_tmp/N;
                 
       end
   end
   
   rho_1 = (MASS_HE*c_He_1 + MASS_O*c_O_1 + MASS_N2*c_N2_1);
   rho_2 = (MASS_HE*c_He_2 + MASS_O*c_O_2 + MASS_N2*c_N2_2); 
   
   subplot(2,2,1)
   hold on
   semilogy(180*lat/pi, c_N2_2./c_N2_1);
   
   subplot(2,2,2)
   hold on
   semilogy(180*lat/pi, c_O_2./c_O_1);
   
   subplot(2,2,3)
   hold on
   semilogy(180*lat/pi, rho_2./rho_1);
   
   subplot(2,2,4)
   hold on
   semilogy(180*lat/pi, c_He_2./c_He_1);
   
   z = 600;
   
   c_He_1 = zeros(1,N);
   c_O_1 = zeros(1,N);
   c_N2_1 = zeros(1,N);

   c_He_2 = zeros(1,N);
   c_O_2 = zeros(1,N);
   c_N2_2 = zeros(1,N);
   
   for i = 1:N
       [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
       
       for j = 1:N
           [~, ~, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
               z, coeffs77, P, f_1, f_1, k_1, t(j), d);
           
           c_He_1(i) = c_He_1(i) + c_He_tmp/N;
           c_O_1(i) = c_O_1(i) + c_O_tmp/N;
           c_N2_1(i) = c_N2_1(i) + c_N2_tmp/N;
           
           [~, ~, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
               z, coeffs77, P, f_2, f_2, k_2, t(j), d);
         
          
           c_He_2(i) = c_He_2(i) + c_He_tmp/N;
           c_O_2(i) = c_O_2(i) + c_O_tmp/N;
           c_N2_2(i) = c_N2_2(i) + c_N2_tmp/N;
                 
       end
   end
   
   rho_1 = (MASS_HE*c_He_1 + MASS_O*c_O_1 + MASS_N2*c_N2_1);
   rho_2 = (MASS_HE*c_He_2 + MASS_O*c_O_2 + MASS_N2*c_N2_2); 
   
   subplot(2,2,1)
   hold on
   semilogy(180*lat/pi, c_N2_2./c_N2_1);
   
   subplot(2,2,2)
   hold on
   semilogy(180*lat/pi, c_O_2./c_O_1);
   
   subplot(2,2,3)
   hold on
   semilogy(180*lat/pi, rho_2./rho_1);
   
   subplot(2,2,4)
   hold on
   semilogy(180*lat/pi, c_He_2./c_He_1);
   
   z = 800;
   
   c_He_1 = zeros(1,N);
   c_O_1 = zeros(1,N);
   c_N2_1 = zeros(1,N);

   c_He_2 = zeros(1,N);
   c_O_2 = zeros(1,N);
   c_N2_2 = zeros(1,N);
   
   for i = 1:N
       [P, ~] = associated_legendre(6, sin(lat(i)), 'positive');
       
       for j = 1:N
           [~, ~, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
               z, coeffs77, P, f_1, f_1, k_1, t(j), d);
           
           c_He_1(i) = c_He_1(i) + c_He_tmp/N;
           c_O_1(i) = c_O_1(i) + c_O_tmp/N;
           c_N2_1(i) = c_N2_1(i) + c_N2_tmp/N;
           
           [~, ~, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
               z, coeffs77, P, f_2, f_2, k_2, t(j), d);
         
          
           c_He_2(i) = c_He_2(i) + c_He_tmp/N;
           c_O_2(i) = c_O_2(i) + c_O_tmp/N;
           c_N2_2(i) = c_N2_2(i) + c_N2_tmp/N;
                 
       end
   end
   
   rho_1 = (MASS_HE*c_He_1 + MASS_O*c_O_1 + MASS_N2*c_N2_1);
   rho_2 = (MASS_HE*c_He_2 + MASS_O*c_O_2 + MASS_N2*c_N2_2); 
   
   subplot(2,2,1)
   hold on
   semilogy(180*lat/pi, c_N2_2./c_N2_1);
   ylim([0.1, 110]);
   xlim([0,90]);
   xlabel('Latitude');
   ylabel('Ratio "disturbed"/"quiet"');
   legend('200 km', '400 km', '600 km', '800 km');
   title('N_2 concentration (particles/cm^3)', 'interpreter', 'latex');
   
   subplot(2,2,2)
   hold on
   semilogy(180*lat/pi, c_O_2./c_O_1);
   ylim([0.1, 110]);
   xlim([0,90]);
   xlabel('Latitude');
   ylabel('Ratio "disturbed"/"quiet"');
   title('O concentration (particles/cm^3)', 'interpreter', 'latex');
   
   subplot(2,2,3)
   hold on
   semilogy(180*lat/pi, rho_2./rho_1);
   ylim([0.1, 10]);
   xlim([0,90]);
   xlabel('Latitude');
   ylabel('Ratio "disturbed"/"quiet"');
   title('Density (g/cm^3)', 'interpreter', 'latex');
   
   subplot(2,2,4)
   hold on
   semilogy(180*lat/pi, c_He_2./c_He_1);
   ylim([0.1, 10]);
   xlim([0,90]);
   xlabel('Latitude');
   ylabel('Ratio "disturbed"/"quiet"');
   title('He concentration (particles/cm^3)', 'interpreter', 'latex');
   
end

%% figure 6
% Sensitivity of density and concentrations to latitude, solar time and day
% number during high solar activty
if max(FIGS_TO_PLOT == 6)

    N = 100;
    M = 100;
    lat = linspace(-pi/2, pi/2, N);
    
    t = 15;
    d =linspace(0, 365, M);
    
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
% Sensitivity of density and concentrations to latitude, solar time and day
% number during low solar activty
if max(FIGS_TO_PLOT == 7)
    
    N = 100;
    M = 100;
    lat = linspace(-pi/2, pi/2, N);
    
    t = 15;
    d =linspace(0, 365, M);
    
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
% Sensitivity of density and concentrations to solar time
if max(FIGS_TO_PLOT == 8)
    
    N = 30;
    d = linspace(0,365,N);
    k = 1;
    lat = 0;
    z = 400;
    f = 92;   
    
    t = linspace(0,24,N);
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    
    T = zeros(1,N);
    c_He = zeros(1,N);
    c_O = zeros(1,N);
    c_N2 = zeros(1,N);
    
    for i = 1:N        
        for j = 1:N
            
            [~, T_tmp, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
                    z, coeffs77, P, f, f, k, t(i), d(j));
                
            T(i) = T(i) + T_tmp/N;
            c_He(i) = c_He(i) + c_He_tmp/N;
            c_O(i) = c_O(i) + c_O_tmp/N;
            c_N2(i) = c_N2(i) + c_N2_tmp/N;    
        end
    end
    
    rho = (MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    figure();
    subplot(4,2,2);
    plot(t, c_N2);
    xlabel('Time');
    ylim([0, 8e6]);
    xlim([0, 24]);
    annotation('textbox', [0 0.89 1 0.1], 'String', 'N_2 concentration (particles/cm^3)', ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    subplot(4,2,4);
    plot(t, c_O);
    xlabel('Time');
    ylim([0,8e7]);
    xlim([0, 24]);
    annotation('textbox', [0 0.671 1 0.1], 'String', 'O concentration (particles/cm^3)',...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    subplot(4,2,6);
    plot(t, c_He);
    xlabel('Time');
    ylim([0,6e6]);
    xlim([0, 24]);
    annotation('textbox', [0 0.45 1 0.1], 'String', 'He concentration (particles/cm^3)',...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
    
    subplot(4,2,8);
    plot(t, rho);
    xlabel('Time');
    ylim([0,3e-15]);
    xlim([0, 24]);
    annotation('textbox', [0 0.225 1 0.1], 'String', 'Total density (g/cm^3)', ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Fontsize', 15);
 
    k = 2;
    f = 150;
    
    T = zeros(1,N);
    c_He = zeros(1,N);
    c_O = zeros(1,N);
    c_N2 = zeros(1,N);
    
    for i = 1:N
        for j = 1:N
            
            [~, T_tmp, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
                    z, coeffs77, P, f, f, k, t(i), d(j));
                
            T(i) = T(i) + T_tmp/N;
            c_He(i) = c_He(i) + c_He_tmp/N;
            c_O(i) = c_O(i) + c_O_tmp/N;
            c_N2(i) = c_N2(i) + c_N2_tmp/N;    
        end
    end
    
    rho = (MASS_HE*c_He + MASS_O*c_O + MASS_N2*c_N2);
    
    subplot(4,2,1);
    plot(t, c_N2);
    xlabel('Time');
    ylim([0, 2e7]);
    xlim([0, 24]);
    
    subplot(4,2,3);
    plot(t, c_O);
    xlabel('Time');
    ylim([0,3e8]);
    xlim([0, 24]);
    
    subplot(4,2,5);
    plot(t, c_He);
    xlabel('Time');
    ylim([0,8e6]);
    xlim([0, 24]);
    
    subplot(4,2,7);
    plot(t, rho);
    xlabel('Time');
    ylim([0,8e-15]);  
    xlim([0, 24]);
end

%% figures 9, 10, 11 and 12
% Figures 9, 10, 11 and 12 use the same data, varying day number
% Figure 9 is sensitivity of density 
% Figure 10 is sensitivity of thermopause temperature and N_2
% concentration
% Figure 11 is sensitivity of O
% Figure 12 is sensitivity of He
if (max(FIGS_TO_PLOT == 9) || max(FIGS_TO_PLOT == 10) ||...
        max(FIGS_TO_PLOT == 11) || max(FIGS_TO_PLOT == 12))
    
    N = 60;
    t = linspace(0,24,N);
    f = 150;
    k = 2;
    lat = 45*pi/180;
    z = 400;    
    
    d = linspace(0, 365, N);
    
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    
    T_inf_p45 = zeros(1,N);
    c_He_p45 = zeros(1,N);
    c_O_p45 = zeros(1,N);
    c_N2_p45 = zeros(1,N);
    
    for i = 1:N
        for j = 1:N
            [T_tmp, ~, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
                    z, coeffs77, P, f, f, k, t(j), d(i));
                
            T_inf_p45(i) = T_inf_p45(i) + T_tmp/N;
            c_He_p45(i) = c_He_p45(i) + c_He_tmp/N;
            c_O_p45(i) = c_O_p45(i) + c_O_tmp/N;
            c_N2_p45(i) = c_N2_p45(i) + c_N2_tmp/N;   
        end
    end
    
    rho_p45 = (MASS_HE*c_He_p45 + MASS_O*c_O_p45 + MASS_N2*c_N2_p45);
    
    lat = -45*pi/180;
    [P, ~] = associated_legendre(6, sin(lat), 'positive');
    
    T_inf_m45 = zeros(1,N);
    c_He_m45 = zeros(1,N);
    c_O_m45 = zeros(1,N);
    c_N2_m45 = zeros(1,N);
    
    for i = 1:N
              for j = 1:N
            [T_tmp, ~, c_He_tmp, c_O_tmp, c_N2_tmp] = dtm77(...
                    z, coeffs77, P, f, f, k, t(j), d(i));
                
            T_inf_m45(i) = T_inf_m45(i) + T_tmp/N;
            c_He_m45(i) = c_He_m45(i) + c_He_tmp/N;
            c_O_m45(i) = c_O_m45(i) + c_O_tmp/N;
            c_N2_m45(i) = c_N2_m45(i) + c_N2_tmp/N;   
        end
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
    legend('45^o latitude', '-45^o latitude');
    
    figure();
    hold on
    plot(d, T_inf_p45);
    plot(d, T_inf_m45, 'r');
    
    legend('45^o latitude', '-45^o latitude');
    xlabel('day');
    
    title('Thermopause temperature (^oK)', 'Fontsize', 15);
    xlim([min(d), max(d)]);
    
end

%% figure 13
% Sensitivity of He concentration to latitude
if max(FIGS_TO_PLOT == 13)
    
    t = 15; 
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

