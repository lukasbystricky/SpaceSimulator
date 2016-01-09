addpath ..
addpath ../coefficients
addpath ../../../Utilities/Math

close all

header = [  3   171   231   309   342   366   387   405   423   441   753   819   899;
    14    10    13    11     8     7     6     6     6    13    11    10    10;
     4     2     2     1     1     1     1     1     1     8     2     4     4];
 
ephemeris = read_ephemeris_file('ephemeris1960-2200.txt', header);

t = linspace(ephemeris.time_intervals(1), ephemeris.time_intervals(end), 10000);

posM = calculate_postitions(ephemeris, t, 'mercury');
posV = calculate_postitions(ephemeris, t, 'venus');
posE = calculate_postitions(ephemeris, t, 'earth');
posMa = calculate_postitions(ephemeris, t, 'mars');
posJ = calculate_postitions(ephemeris, t, 'jupiter');
posS = calculate_postitions(ephemeris, t, 'saturn');
posU = calculate_postitions(ephemeris, t, 'uranus');
posN = calculate_postitions(ephemeris, t, 'neptune');
posP = calculate_postitions(ephemeris, t, 'pluto');
posSun = calculate_postitions(ephemeris, t, 'sun');
posMoon = calculate_postitions(ephemeris, t, 'moon');

massM = 3.285e23;
massV = 4.867e24;
massMa = 6.4171e23;
massJ = 1.89813e27;
massS = 5.683e26;
massU = 8.681e25;
massN = 1.024e26;
massP = 1.30900e22;
massSun = 1.989e30;
massMoon = 7.347673e22;

figure();
ha = gca;

figure();
hd = gca;

[d, a] = calculate_acceleration(posSun(1:80,:), posE(1:80,:), massSun, false);
semilogy(ha, t(1:80), a, 'linewidth', 2);
hold(ha,'on');
semilogy(hd, t(1:80), d, 'linewidth', 2);
hold(hd,'on');

[d, a] = calculate_acceleration(posMoon(1:80,:), posE(1:80,:), massMoon, true);
semilogy(ha, t(1:80), a, 'linewidth', 2);
semilogy(hd, t(1:80), d, 'linewidth', 2);

[d, a] = calculate_acceleration(posM(1:80,:), posE(1:80,:), massM, false);
semilogy(ha, t(1:80), a, 'linewidth', 2);
semilogy(hd, t(1:80), d, 'linewidth', 2);

[d, a] = calculate_acceleration(posV(1:80,:), posE(1:80,:), massV, false);
semilogy(ha, t(1:80), a, 'linewidth', 2);
semilogy(hd, t(1:80), d, 'linewidth', 2);

[d, a] = calculate_acceleration(posMa(1:80,:), posE(1:80,:), massMa, false);
semilogy(ha, t(1:80), a, 'linewidth', 2);
semilogy(hd, t(1:80), d, 'linewidth', 2);

ha.XTick = [t(10), t(70)];
ha.XTickLabel = {'26 Feb, 1960', '08 Aug, 1961'};

hd.XTick = [t(10), t(70)];
hd.XTickLabel = {'26 Feb, 1960', '08 Aug, 1961'};

xlim(ha, [t(1), t(80)]);
xlim(hd, [t(1), t(80)]);

title(ha, 'Inner bodies gravitational acceleration on Earth [m/s^2]');
title(hd, 'Inner bodies distance from Earth [m]');

legend(ha,'Sun', 'Moon', 'Mercury', 'Venus', 'Mars');
legend(hd,'Sun', 'Moon', 'Mercury', 'Venus', 'Mars');

figure();
ha = gca;

figure();
hd = gca;

[d, a] = calculate_acceleration(posJ, posE, massJ, false);
semilogy(ha, t, a);
hold(ha, 'on');
semilogy(hd, t, d);
hold(hd, 'on');

[d, a] = calculate_acceleration(posS, posE, massS, false);
semilogy(ha, t, a);
semilogy(hd, t, d);

[d, a] = calculate_acceleration(posU, posE, massU, false);
semilogy(ha, t, a);
semilogy(hd, t, d);

[d, a] = calculate_acceleration(posN, posE, massN, false);
semilogy(ha, t, a);
semilogy(hd, t, d);

[d, a] = calculate_acceleration(posP, posE, massP, false);

semilogy(ha, t, a);
semilogy(hd, t, d);

ha.XTick = [t(1000), t(end - 1000)];
ha.XTickLabel = {'14 Jan, 1984', '5 Jan, 2177'};

hd.XTick = [t(1000), t(end - 1000)];
hd.XTickLabel = {'14 Jan, 1984', '5 Jan, 2117'};

xlim(ha, [t(1), t(end)]);
xlim(hd, [t(1), t(end)]);

title(ha, 'Outer planets gravitational acceleration on Earth [m/s^2]');
title(hd, 'Outer planets distance from Earth [m]');

legend(ha,'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto');
legend(hd,'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto');