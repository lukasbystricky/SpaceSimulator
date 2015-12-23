az = linspace(0,2*pi,80);
el = linspace(-pi/2,pi/2,80);

[AZ, EL] = meshgrid(az, el);

n = 20;
m = 20;

u = zeros(size(AZ));

for i = 1:length(az)
    for j = 1:length(el)
        u(i,j) = spheric_harmonic(AZ(i,j), EL(i,j), n, m);
    end
end

figure();
[X, Y, Z] = sph2cart(AZ, EL, r);
surf(X, Y, Z, u, 'EdgeColor', 'none');
shading interp
colorbar;
