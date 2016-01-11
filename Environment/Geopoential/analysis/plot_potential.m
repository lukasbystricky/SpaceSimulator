addpath ../
addpath ../coefficients
addpath ../../../Utilities/Math

close all

MU = 3.986004415e14;
R_EQ = 6.3781363e6;

N = 100;
M_MAX = 5;
az = linspace(-pi, pi,N);
el = linspace(-pi/2,pi/2, N);

[AZ, EL] = meshgrid(az, el(end:-1:1));
r = R_EQ + 200e3;                        %200 km above mean Earth surface

u = zeros(size(AZ));
g = zeros(size(AZ));

coeffs_all = dlmread('egm96_to360.ascii');
m_rows = find(coeffs_all(:,1) == M_MAX);
coeffs = coeffs_all(1:m_rows(1),:);

for i = 1:length(az)
    for j = 1:length(el)
        
        [u(i,j), g_vec] = geopotential(r, AZ(i,j), EL(i,j), coeffs);
        
        
        
        g_spherical = [-MU/r^2; 0; 0];
        
        rot = [sin(EL(i,j))*cos(AZ(i,j)), cos(EL(i,j))*cos(AZ(i,j)), -sin(AZ(i,j));
                sin(EL(i,j))*sin(AZ(i,j)), cos(EL(i,j))*sin(AZ(i,j)), cos(AZ(i,j));
                cos(EL(i,j)),         -sin(EL(i,j)),        0];

        g_spherical = rot*g_spherical;
        
        g(i,j) = norm(g_vec - g_spherical);
    end
end

%% surface plot on sphere of potential and acceleration
figure();
[X, Y, Z] = sph2cart(AZ, EL, r);
surf(X, Y, Z, u, 'EdgeColor', 'none');
title('Gravitational Potential');
colorbar;
shading interp

figure();
surf(X, Y, Z, g, 'EdgeColor', 'none');
title('Gravitational Acceleration');
colorbar;
shading interp


%% map of earth acceleration
figure();

img = imread('coast.jpg');
imagesc([-180,180],[-90, 90], flipud(img));
hold all;

[~,h] = contourf(az*(180/pi), el*(180/pi), flipud(g), 30, 'EdgeColor', 'none');

set(gca, 'clim',[min(min(g)), max(max(g))]);
set(gca,'YDir','normal');

colormap(jet)
colorbar

pause(0.05)
hFaces = h.FacePrims;
for faceIdx = 1 : numel(hFaces)
   hFaces(faceIdx).ColorType = 'truecoloralpha';  % default = 'truecolor'
   hFaces(faceIdx).ColorData(4) = 150;   % default=255
end

title('Perturbing Acceleration [m/s^2]');
