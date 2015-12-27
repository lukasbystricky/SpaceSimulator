addpath ../
addpath ../coefficients
addpath ../../../Utilities/Math

close all

N = 100;
az = linspace(0,2*pi,N);
el = linspace(-pi/2,pi/2,N);

[AZ, EL] = meshgrid(az, el);
r = 6.371e6 + 200e3*1.6; %200 miles above mean Earth surface
mu = 3.986004415e14;

u = zeros(size(AZ));
g = zeros(size(AZ));

coeffs_all = dlmread('egm96_to360.ascii');
coeffs = coeffs_all(2:50,:);

for i = 1:length(az)
    for j = 1:length(el)
        [u(i,j), g_vec] = geopotential(r, AZ(i,j), EL(i,j), coeffs);
        
        g(i,j) = norm(g_vec);        
    end
end

figure();
[X, Y, Z] = sph2cart(AZ, EL, r);
surf(X, Y, Z, u, 'EdgeColor', 'none');
title('Gravitational Potential');
colorbar;
shading interp

figure;
surf(X, Y, Z, g - mu/r^2, 'EdgeColor', 'none');
title('Gravitational Acceleration');
colorbar;
shading interp

figure;

gtmp1 = g(:,1:N/2);
gtmp2 = g(:,N/2+1:end);

gfinal = [gtmp2, gtmp1];

img = imread('coast.jpg');
imagesc([-180,180],[-90, 90], img);
hold all;

[~,h] = contourf(az*(180/pi)-180, el*(180/pi), flipud(gfinal- mu/r^2), 30, 'EdgeColor', 'none');
set(gca, 'clim',[min(min(gfinal- mu/r^2)), max(max(gfinal- mu/r^2))]);

colormap(jet)
colorbar

pause(0.05)
hFaces = h.FacePrims;
for faceIdx = 1 : numel(hFaces)
   hFaces(faceIdx).ColorType = 'truecoloralpha';  % default = 'truecolor'
   hFaces(faceIdx).ColorData(4) = 150;   % default=255
end