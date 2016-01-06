addpath ../
addpath ../coefficients
addpath ../../../Utilities/Math

close all

N = 100;
M_MAX = 15;
az = linspace(0,2*pi,N);
el = linspace(-pi/2,pi/2,N);

[AZ, EL] = meshgrid(az, el);

MU = 3.986004415e14;
%R_EQ = 6.378137e6;
R_EQ = 6.3781363e6;
R_POL = 6.356752314245e6;

n = zeros(size(AZ));
t = zeros(size(AZ));

coeffs_all = dlmread('egm96_to360.ascii');

m_rows = find(coeffs_all(:,1) == M_MAX);
coeffs = coeffs_all(1:m_rows(1),:);

zonal_rows = find(coeffs(:,2) == 0);
coeffs_zonal = coeffs(zonal_rows, :);

for i = 1:length(az)
    for j = 1:length(el)
        
        %r = sqrt(((R_EQ^2*cos(EL(i,j)))^2 + (R_POL^2*sin(EL(i,j)))^2)/((R_EQ*cos(EL(i,j)))^2 + (R_POL*sin(EL(i,j)))^2));
        r = R_EQ;
        [u, ~] = geopotential(r, AZ(i,j), EL(i,j), coeffs);
        
        [~, n(i,j)] = calculate_undulation(r, AZ(i,j), EL(i,j), u, coeffs_zonal, R_EQ, MU);    
    end
end

figure();

n_tmp1 = n(:,1:N/2);
n_tmp2 = n(:,N/2+1:end);

n_final = [n_tmp2, n_tmp1];

img = imread('coast.jpg');
imagesc([-180,180],[-90, 90], img);
hold all;

[~,h] = contourf(az*(180/pi)-180, el*(180/pi), -flipud(n_final), 30, 'EdgeColor', 'none');
%set(gca, 'clim',[min(min(n_final)), max(max(n_final))]);
set(gca, 'clim',[-110, 85]);
shading interp;

colormap(jet)
colorbar

pause(0.05)
hFaces = h.FacePrims;
for faceIdx = 1 : numel(hFaces)
   hFaces(faceIdx).ColorType = 'truecoloralpha';  % default = 'truecolor'
   hFaces(faceIdx).ColorData(4) = 150;   % default=255
end

title('Geoid Undulation');
