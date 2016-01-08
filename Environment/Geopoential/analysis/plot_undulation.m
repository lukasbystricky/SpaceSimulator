%% This script does overlays a contour plot of the Earth's undulation on a
% map of the Earth

addpath ../
addpath ../coefficients
addpath ../../../Utilities/Math

close all

N = 100;
M_MAX = 15;
az = linspace(-pi,pi,N);
el = linspace(-pi/2,pi/2,N);

[AZ, EL] = meshgrid(az, el(end:-1:1));

MU = 3.986004415e14;
R_EQ = 6.3781363e6;

n = zeros(size(AZ));
t = zeros(size(AZ));

coeffs_all = dlmread('egm96_to360.ascii');

m_rows = find(coeffs_all(:,1) == M_MAX);
coeffs = coeffs_all(1:m_rows(1),:);

zonal_rows = find(coeffs(:,2) == 0);
coeffs_zonal_ = coeffs(zonal_rows, :);

for i = 1:length(az)
    for j = 1:length(el)
        
        [u, ~] = geopotential(R_EQ, AZ(i,j), EL(i,j), coeffs);
        
        [~, n(i,j)] = calculate_undulation(R_EQ, AZ(i,j), EL(i,j), u, coeffs_zonal(1:2:end,:), R_EQ, MU);    
    end
end

figure();

img = imread('coast.jpg');
image([-180,180],[-90, 90], flipud(img));

hold all;
[~,h] = contourf(az*(180/pi), el*(180/pi), -flipud(n), 30, 'EdgeColor', 'none');

shading interp;

set(gca, 'clim',[-110, 85]);
set(gca,'YDir','normal');

colormap(jet)
colorbar

pause(0.05)
hFaces = h.FacePrims;
for faceIdx = 1 : numel(hFaces)
   hFaces(faceIdx).ColorType = 'truecoloralpha';  % default = 'truecolor'
   hFaces(faceIdx).ColorData(4) = 150;   % default=255
end

title('Geoid Undulation [m]');
