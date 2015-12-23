N = 50;
az = linspace(0,2*pi,N);
el = linspace(-pi/2,pi/2,N);

[AZ, EL] = meshgrid(az, el);
r = 6.371e6 + 200e3*1.6; %200 miles above mean Earth surface
mu = 3.986004415e14;

% J = [0, 1.7555e25, 0];
% %J = [0, 0, -2.619e29];
% S = 0;
% C = 0;

u = zeros(size(AZ));
g = zeros(size(AZ));

coeffs_all = dlmread('egm96_to360.ascii');
coeffs = coeffs_all(2:300,:);

for i = 1:length(az)
    for j = 1:length(el)
        [u(i,j), g_vec] = gravitational_acceleration_harmonics(r, AZ(i,j), EL(i,j), coeffs);
        
        g(i,j) = norm(g_vec);        
    end
end

figure();
[X, Y, Z] = sph2cart(AZ, EL, r);
surf(X, Y, Z, u, 'EdgeColor', 'none');
colorbar;
shading interp

figure;
surf(X, Y, Z, g - mu/r^2, 'EdgeColor', 'none');
colorbar;
shading interp

opengl software;
figure;

gtmp1 = g(:,1:N/2);
gtmp2 = g(:,N/2+1:end);

gfinal = [gtmp2, gtmp1];

img = imread('coast.jpg');
imagesc([-180,180],[-90, 90], img);
hold all;

[C,h] = contourf(az*(180/pi)-180, el*(180/pi), flipud(gfinal- mu/r^2), 30, 'EdgeColor', 'none');
set(gca, 'clim',[min(min(gfinal- mu/r^2)), max(max(gfinal- mu/r^2))]);
%set(gcf, 'Renderer', 'OpenGL');
colorbar

ph = findobj(h, '-property', 'FaceAlpha');
set(ph, 'FaceAlpha', 0.05);
% for k = [1:-0.1:0.1 0.1:0.1:1]
%  set(ph, 'FaceAlpha', k);
%  pause(0.5)
% end






