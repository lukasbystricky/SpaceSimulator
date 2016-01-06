addpath ../
addpath ../coefficients
addpath ../../../Utilities/Math

close all

N = 100;
az = linspace(0,2*pi,N);
el = linspace(-pi/2,pi/2,N);

[AZ, EL] = meshgrid(az, el);
%r = 6.371e6 + 200e3*1.6; %200 miles above mean Earth surface
MU = 3.986004415e14;
R_EQ = 6.378137e6;
R_POL = 6.356752314245e6;

r_expected = zeros(size(AZ));
u = zeros(size(AZ));
g = zeros(size(AZ));

coeffs_all = dlmread('egm96_to360.ascii');
coeffs = coeffs_all(2:100,:);

for i = 1:length(az)
    for j = 1:length(el)
        
        r_expected(i,j) = sqrt(((R_EQ^2*cos(EL(i,j)))^2 + (R_POL^2*sin(EL(i,j)))^2)/((R_EQ*cos(EL(i,j)))^2 + (R_POL*sin(EL(i,j)))^2));
        [u(i,j), g_vec] = geopotential(r_expected(i,j), AZ(i,j), EL(i,j), coeffs);
        
        g(i,j) = norm(g_vec);        
    end
end

%% surface plot on sphere of potential and acceleration
figure();
[X, Y, Z] = sph2cart(AZ, EL, r);
surf(X, Y, Z, u, 'EdgeColor', 'none');
title('Gravitational Potential');
colorbar;
shading interp

figure;
surf(X, Y, Z, g - MU/r^2, 'EdgeColor', 'none');
title('Gravitational Acceleration');
colorbar;
shading interp


%% map of earth acceleration
figure();

gtmp1 = g(:,1:N/2);
gtmp2 = g(:,N/2+1:end);

gfinal = [gtmp2, gtmp1];

img = imread('coast.jpg');
imagesc([-180,180],[-90, 90], img);
hold all;

[~,h] = contourf(az*(180/pi)-180, el*(180/pi), flipud(gfinal- MU/r^2), 30, 'EdgeColor', 'none');
set(gca, 'clim',[min(min(gfinal- MU/r^2)), max(max(gfinal- MU/r^2))]);

colormap(jet)
colorbar

pause(0.05)
hFaces = h.FacePrims;
for faceIdx = 1 : numel(hFaces)
   hFaces(faceIdx).ColorType = 'truecoloralpha';  % default = 'truecolor'
   hFaces(faceIdx).ColorData(4) = 150;   % default=255
end

title('Gravitational Acceleration');

%% map of height
figure();

r_actual = zeros(size(gfinal));

height_anomaly = zeros(size(gfinal));

for i = 1:length(el)
    for j = 1:length(az)
        %r_actual(i,j) = sqrt(MU/gfinal(i,j));
        height_anomaly(i,j) = (u(i,j) + MU/r_expected(i,j))/r_expected(i,j);
    end
end

imagesc([-180,180],[-90, 90], img);
hold all;

% [~,h] = contourf(az*(180/pi)-180, el*(180/pi), -flipud(r_actual - r_expected), 30, 'EdgeColor', 'none');
% set(gca, 'clim',[min(min(r_actual - r_expected)), max(max(r_actual - r_expected))]);

[~,h] = contourf(az*(180/pi)-180, el*(180/pi), -flipud(height_anomaly), 30, 'EdgeColor', 'none');
set(gca, 'clim',[min(min(height_anomaly)), max(max(height_anomaly))]);

colormap(jet)
colorbar

pause(0.05)
hFaces = h.FacePrims;
for faceIdx = 1 : numel(hFaces)
   hFaces(faceIdx).ColorType = 'truecoloralpha';  % default = 'truecolor'
   hFaces(faceIdx).ColorData(4) = 150;   % default=255
end

title('Height');


