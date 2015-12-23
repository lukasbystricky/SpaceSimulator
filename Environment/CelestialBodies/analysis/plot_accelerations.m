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

massM = 3.285e23;
massV = 4.867e24;
massMa = 6.39e23;
massJ = 1.89813e27;
massS = 5.683e26;
massU = 8.681e25;
massN = 1.024e26;
massP = 1.30900e22;
massSun = 1.989e30;


figure();
ha = gca;

figure();
hd = gca;


[d, a] = calculate_acceleration(posSun, posE, massSun);
semilogy(ha, a);
hold(ha,'on')
semilogy(hd, d);
hold(hd,'on')

[d, a] = calculate_acceleration(posM, posE, massM);
semilogy(ha, a);
semilogy(hd, d);

[d, a] = calculate_acceleration(posV, posE, massV);
semilogy(ha, a);
semilogy(hd, d);

[d, a] = calculate_acceleration(posMa, posE, massMa);
semilogy(ha, a);
semilogy(hd, d);

[d, a] = calculate_acceleration(posJ, posE, massJ);
semilogy(ha, a);
semilogy(hd, d);

[d, a] = calculate_acceleration(posS, posE, massS);
semilogy(ha, a);
semilogy(hd, d);

[d, a] = calculate_acceleration(posU, posE, massU);
semilogy(ha, a);
semilogy(hd, d);

[d, a] = calculate_acceleration(posN, posE, massN);
semilogy(ha, a);
semilogy(hd, d);

[d, a] = calculate_acceleration(posP, posE, massP);
semilogy(ha, a);
semilogy(hd, d);

legend(ha,'sun', 'mercury', 'venus', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune', 'pluto');
legend(hd,'sun', 'mercury', 'venus', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune', 'pluto');