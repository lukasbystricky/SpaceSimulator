%% Create an animated gif of the orbits of the inner planets

addpath ..
addpath ../coefficients
addpath ../../../Utilities/Math

header = [  3   171   231   309   342   366   387   405   423   441   753   819   899;
    14    10    13    11     8     7     6     6     6    13    11    10    10;
     4     2     2     1     1     1     1     1     1     8     2     4     4];


ephemeris = read_ephemeris_file('ephemeris1940-2220.txt', header);
t = linspace(ephemeris.time_intervals(1), ephemeris.time_intervals(25), 1000);

posMerc = calculate_postitions(ephemeris, t, 'mercury');
posVen = calculate_postitions(ephemeris, t, 'venus');
posEarth = calculate_postitions(ephemeris, t, 'earth');
posMars = calculate_postitions(ephemeris, t, 'mars');
posSun = calculate_postitions(ephemeris, t, 'sun');

h = figure();

view(10,20);
xmax = 1.15*max(abs(posMars(:,1)));
ymax = 1.15*max(abs(posMars(:,1)));
zmax = 1.15*max(abs(posMars(:,1)));

for i = 1:10:length(t)
   
   hold off
   plot3(posSun(i,1), posSun(i,2), posSun(i,3), 'yo');
   
   xlim([-xmax, xmax]);
   ylim([-ymax, ymax]);
   zlim([-zmax, zmax]);
   
   set(gca, 'color', [0 0 0]);
   axis off
   set(gcf, 'color', [0 0 0]);
   
   hold on;
   
   plot3(posMerc(i,1), posMerc(i,2), posMerc(i,3), 'bo');
   plot3(posVen(i,1), posVen(i,2), posVen(i,3), 'ro');
   plot3(posEarth(i,1), posEarth(i,2), posEarth(i,3), 'go');
   plot3(posMars(i,1), posMars(i,2), posMars(i,3), 'mo');
   
   
   legend({'Sun', 'Mercury','Venus','Earth','Mars'}, 'TextColor', 'w', 'location', 'northwest');
   
   plot3(posMerc(1:i,1), posMerc(1:i,2), posMerc(1:i,3), '-b', 'linewidth', 2);
   plot3(posVen(1:i,1), posVen(1:i,2), posVen(1:i,3), '-r', 'linewidth', 2);
   plot3(posEarth(1:i,1), posEarth(1:i,2), posEarth(1:i,3), '-g', 'linewidth', 2);
   plot3(posMars(1:i,1), posMars(1:i,2), posMars(1:i,3), '-m', 'linewidth', 2);
   plot3(posSun(1:i,1), posSun(1:i,2), posSun(1:i,3), '-y', 'linewidth', 2);   
   
   view([37, 90]);
   
   frame = getframe;
   im = frame2im(frame);
   [A,map] = rgb2ind(im,256); 

   
   if (i == 1)       
       imwrite(A,map,'orbit_inner.gif','gif','LoopCount',Inf,'DelayTime',1);
   else       
       imwrite(A,map,'orbit_inner.gif','gif','WriteMode','append','DelayTime', 0);
   end
end