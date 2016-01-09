%% Create an animated gif of the orbits of the outer planets

addpath ..
addpath ../coefficients
addpath ../../../Utilities/Math

header = [  3   171   231   309   342   366   387   405   423   441   753   819   899;
    14    10    13    11     8     7     6     6     6    13    11    10    10;
     4     2     2     1     1     1     1     1     1     8     2     4     4];
 
ephemeris = read_ephemeris_file('ephemeris1940-2220.txt', header);
t = linspace(ephemeris.time_intervals(1), ephemeris.time_intervals(end), 1000);

posJ = calculate_postitions(ephemeris, t, 'jupiter');
posS = calculate_postitions(ephemeris, t, 'saturn');
posU = calculate_postitions(ephemeris, t, 'uranus');
posN = calculate_postitions(ephemeris, t, 'neptune');
posP = calculate_postitions(ephemeris, t, 'pluto');

h = figure();

view(10,20);
xmax = 1.05*max([abs(posN(:,1)); abs(posP(:,1))]);
ymax = 1.05*max([abs(posN(:,2)); abs(posP(:,2))]);
zmax = 1.05*max([abs(posN(:,3)); abs(posP(:,3))]);

for i = 1:10:length(t);
   
   hold off
   plot3(posJ(i,1), posJ(i,2), posJ(i,3), 'bo');
   xlim([-xmax, xmax]);
   ylim([-ymax, ymax]);
   zlim([-zmax, zmax]);
   
   set(gca, 'color', [0 0 0]);
   axis off
   set(gcf, 'color', [0 0 0]);
   
   hold on;
   plot3(posS(i,1), posS(i,2), posS(i,3), 'ro');
   plot3(posU(i,1), posU(i,2), posU(i,3), 'go');
   plot3(posN(i,1), posN(i,2), posN(i,3), 'mo');
   plot3(posP(i,1), posP(i,2), posP(i,3), 'yo');
   
   legend({'Jupiter','Saturn','Uranus','Neptune','Pluto'}, 'TextColor', 'w', 'location', 'northwest');
   
   plot3(posJ(1:i,1), posJ(1:i,2), posJ(1:i,3), '-b', 'linewidth', 2);
   plot3(posS(1:i,1), posS(1:i,2), posS(1:i,3), '-r', 'linewidth', 2);
   plot3(posU(1:i,1), posU(1:i,2), posU(1:i,3), '-g', 'linewidth', 2);
   plot3(posN(1:i,1), posN(1:i,2), posN(1:i,3), '-m', 'linewidth', 2);
   plot3(posP(1:i,1), posP(1:i,2), posP(1:i,3), '-y', 'linewidth', 2);   
   
   view([37, 90]);
   
   frame = getframe;
   im = frame2im(frame);
   [A,map] = rgb2ind(im,256); 

   
   if (i == 1)       
       imwrite(A,map,'orbit_outer.gif','gif','LoopCount',Inf,'DelayTime',1);
   else       
       imwrite(A,map,'orbit_outer.gif','gif','WriteMode','append','DelayTime', 0);
   end
end