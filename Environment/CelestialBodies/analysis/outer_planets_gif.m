header = [  3   171   231   309   342   366   387   405   423   441   753   819   899;
    14    10    13    11     8     7     6     6     6    13    11    10    10;
     4     2     2     1     1     1     1     1     1     8     2     4     4];
 
% ephemeris1 = read_ephemeris_file('ascp1960.405', header);
% ephemeris2 = read_ephemeris_file('ascp1980.405', header);
% ephemeris3 = read_ephemeris_file('ascp2000.405', header);

% t1 = linspace(ephemeris1.time_intervals(1), ephemeris1.time_intervals(end), 1000);
% t2 = linspace(ephemeris2.time_intervals(1), ephemeris2.time_intervals(end), 1000);
% t3 = linspace(ephemeris3.time_intervals(1), ephemeris3.time_intervals(end), 1000);

ephemeris1 = read_ephemeris_file('ephemeris1940-2220.txt', header);
t1 = linspace(ephemeris1.time_intervals(1), ephemeris1.time_intervals(end), 10000);

posJ = calculate_postitions(ephemeris1, t1, 'jupiter');
posS = calculate_postitions(ephemeris1, t1, 'saturn');
posU = calculate_postitions(ephemeris1, t1, 'uranus');
posN = calculate_postitions(ephemeris1, t1, 'neptune');
posP = calculate_postitions(ephemeris1, t1, 'pluto');

% posJ2 = calculate_postitions(ephemeris2, t2, 'jupiter');
% posS2 = calculate_postitions(ephemeris2, t2, 'saturn');
% posU2 = calculate_postitions(ephemeris2, t2, 'uranus');
% posN2 = calculate_postitions(ephemeris2, t2, 'neptune');
% posP2 = calculate_postitions(ephemeris2, t2, 'pluto');
% 
% posJ3 = calculate_postitions(ephemeris3, t3, 'jupiter');
% posS3 = calculate_postitions(ephemeris3, t3, 'saturn');
% posU3 = calculate_postitions(ephemeris3, t3, 'uranus');
% posN3 = calculate_postitions(ephemeris3, t3, 'neptune');
% posP3 = calculate_postitions(ephemeris3, t3, 'pluto');
% 
% posJ = [posJ1; posJ2; posJ3];
% posS = [posS1; posS2; posS3];
% posU = [posU1; posU2; posU3];
% posN = [posN1; posN2; posN3];
% posP = [posP1; posP2; posP3];

h = figure();

view(10,20);
xmax = 1.05*max([abs(posN(:,1)); abs(posP(:,1))]);
ymax = 1.05*max([abs(posN(:,2)); abs(posP(:,2))]);
zmax = 1.05*max([abs(posN(:,3)); abs(posP(:,3))]);

for i = 1:10:length(t1);
   
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
   
   %legend('Jupiter','Saturn','Uranus','Neptune','Pluto', 'location', 'south', 'TextColor', 'w' );
   %legend('boxoff');
   
   plot3(posJ(1:i,1), posJ(1:i,2), posJ(1:i,3), '-b', 'linewidth', 2);
   plot3(posS(1:i,1), posS(1:i,2), posS(1:i,3), '-r', 'linewidth', 2);
   plot3(posU(1:i,1), posU(1:i,2), posU(1:i,3), '-g', 'linewidth', 2);
   plot3(posN(1:i,1), posN(1:i,2), posN(1:i,3), '-m', 'linewidth', 2);
   plot3(posP(1:i,1), posP(1:i,2), posP(1:i,3), '-y', 'linewidth', 2);
   
   
   frame = getframe;
   im = frame2im(frame);
   [A,map] = rgb2ind(im,256); 

   
   if (i == 1)       
       imwrite(A,map,'orbit.gif','gif','LoopCount',Inf,'DelayTime',1);
   else       
       imwrite(A,map,'orbit.gif','gif','WriteMode','append','DelayTime', 0);
   end
end