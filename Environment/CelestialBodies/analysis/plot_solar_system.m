
header = [  3   171   231   309   342   366   387   405   423   441   753   819   899;
    14    10    13    11     8     7     6     6     6    13    11    10    10;
     4     2     2     1     1     1     1     1     1     8     2     4     4];
 
ephemeris = read_ephemeris_file('ascp1960.405', header);

t = linspace(ephemeris.time_intervals(1), ephemeris.time_intervals(end), 10000);
figure();
hold

pos = calculate_postitions(ephemeris, t, 'mercury');
plot3(pos(:,1), pos(:,2), pos(:,3));

pos = calculate_postitions(ephemeris, t, 'venus');
plot3(pos(:,1), pos(:,2), pos(:,3));

pos = calculate_postitions(ephemeris, t, 'earth');
plot3(pos(:,1), pos(:,2), pos(:,3));

pos = calculate_postitions(ephemeris, t, 'mars');
plot3(pos(:,1), pos(:,2), pos(:,3));

pos = calculate_postitions(ephemeris, t, 'jupiter');
plot3(pos(:,1), pos(:,2), pos(:,3));

pos = calculate_postitions(ephemeris, t, 'saturn');
plot3(pos(:,1), pos(:,2), pos(:,3));

pos = calculate_postitions(ephemeris, t, 'uranus');
plot3(pos(:,1), pos(:,2), pos(:,3));

pos = calculate_postitions(ephemeris, t, 'neptune');
plot3(pos(:,1), pos(:,2), pos(:,3));

pos = calculate_postitions(ephemeris, t, 'pluto');
plot3(pos(:,1), pos(:,2), pos(:,3));

pos = calculate_postitions(ephemeris, t, 'sun');
plot3(pos(:,1), pos(:,2), pos(:,3));