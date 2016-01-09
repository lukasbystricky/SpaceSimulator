function ephemeris = read_ephemeris_file(fname, header)

%% calculate_positions Reads a planetary ephemeris file
% ephemeris = read_ephemeris_file(fname, header) reads and parses a JPL
% planetary ephemeris file
%
% Given: 
% - a file name
% - a header
%
% this model returns:
% - a ephemeris containing coefficients for each planet/body
%
% The ephemeris file provided must follow the format described at 
% ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/ascii_format.txt. 
% 
% The header matrix must be be exactly the matrix provided in the JPL
% header file. 
%
% The returned ephemeris is a structure with a field for each planet/body
% each containing the a list of the time subintervals and the respective
% coefficients.

M = dlmread(fname);
A = M';

[nrows, ~] = size(M);

n_coeffs = (header(1,end)+header(2,end)*header(3,end)*3 + 1)/3;

n_intervals = floor(nrows/(n_coeffs+1)); % assumes no more than ncoeffs rows

coeffs = zeros(n_intervals,n_coeffs*3 - 2);
times = zeros(n_intervals, 2);

for i = 1:n_intervals    
    times(i,:) = M(i + 1 + (i - 1)*(n_coeffs), 1:2);
    coeffs(i,:) = A(5*i+ (i - 1)*(n_coeffs*3 - 2) + 1: 5*i + i*(n_coeffs*3 - 2));
end

ephemeris.validity = [times(1), times(end)];
ephemeris.time_intervals = times;

for i = 1:11
    [coeffs_x, coeffs_y, coeffs_z, time_subintervals] = parse_coefficients(coeffs, times, header, i);
    
    switch i
        case 1
            
            ephemeris.mercury.time_subintervals = time_subintervals;
            ephemeris.mercury.coeffs_x = coeffs_x;
            ephemeris.mercury.coeffs_y = coeffs_y;
            ephemeris.mercury.coeffs_z = coeffs_z;
            
        case 2
            ephemeris.venus.time_subintervals = time_subintervals;
            ephemeris.venus.coeffs_x = coeffs_x;
            ephemeris.venus.coeffs_y = coeffs_y;
            ephemeris.venus.coeffs_z = coeffs_z;
            
        case 3
            ephemeris.earth.time_subintervals = time_subintervals;
            ephemeris.earth.coeffs_x = coeffs_x;
            ephemeris.earth.coeffs_y = coeffs_y;
            ephemeris.earth.coeffs_z = coeffs_z;
            
        case 4
            ephemeris.mars.time_subintervals = time_subintervals;
            ephemeris.mars.coeffs_x = coeffs_x;
            ephemeris.mars.coeffs_y = coeffs_y;
            ephemeris.mars.coeffs_z = coeffs_z;
            
        case 5
            ephemeris.jupiter.time_subintervals = time_subintervals;
            ephemeris.jupiter.coeffs_x = coeffs_x;
            ephemeris.jupiter.coeffs_y = coeffs_y;
            ephemeris.jupiter.coeffs_z = coeffs_z;
            
        case 6
            ephemeris.saturn.time_subintervals = time_subintervals;
            ephemeris.saturn.coeffs_x = coeffs_x;
            ephemeris.saturn.coeffs_y = coeffs_y;
            ephemeris.saturn.coeffs_z = coeffs_z;
            
        case 7
            ephemeris.uranus.time_subintervals = time_subintervals;
            ephemeris.uranus.coeffs_x = coeffs_x;
            ephemeris.uranus.coeffs_y = coeffs_y;
            ephemeris.uranus.coeffs_z = coeffs_z;
            
        case 8
            ephemeris.neptune.time_subintervals = time_subintervals;
            ephemeris.neptune.coeffs_x = coeffs_x;
            ephemeris.neptune.coeffs_y = coeffs_y;
            ephemeris.neptune.coeffs_z = coeffs_z;
            
        case 9
            ephemeris.pluto.time_subintervals = time_subintervals;
            ephemeris.pluto.coeffs_x = coeffs_x;
            ephemeris.pluto.coeffs_y = coeffs_y;
            ephemeris.pluto.coeffs_z = coeffs_z;
            
        case 10
            ephemeris.moon.time_subintervals = time_subintervals;
            ephemeris.moon.coeffs_x = coeffs_x;
            ephemeris.moon.coeffs_y = coeffs_y;
            ephemeris.moon.coeffs_z = coeffs_z;       
            
        case 11
            ephemeris.sun.time_subintervals = time_subintervals;
            ephemeris.sun.coeffs_x = coeffs_x;
            ephemeris.sun.coeffs_y = coeffs_y;
            ephemeris.sun.coeffs_z = coeffs_z;  
    end
end

end

function [coeffs_x, coeffs_y, coeffs_z, time_intervals] = parse_coefficients(coeffs, times, header, planet_id)

counter = 1;

[nintervals, ~] = size(coeffs);

start = header(1, planet_id);
ncoeffs = header(2, planet_id);
nsubintervals = header(3, planet_id);

time_intervals = zeros(nsubintervals*nintervals,2);
coeffs_x = zeros(nsubintervals*nintervals, ncoeffs);
coeffs_y = zeros(nsubintervals*nintervals, ncoeffs);
coeffs_z = zeros(nsubintervals*nintervals, ncoeffs);

for i = 1:nintervals
    
    coeffs_all = coeffs(i, start - 2 : start - 2 + nsubintervals*ncoeffs*3);
    interval_length = times(i,2) - times(i,1);
    dt = interval_length/nsubintervals;
    
    for j = 1:nsubintervals
        time_intervals(counter,:) = [times(i,1) + (j-1)*dt, times(i,1) + j*dt];
        
        a = 1 + (j-1)*3*ncoeffs;
        coeffs_x(counter,:) = coeffs_all(a : a + ncoeffs - 1);
        coeffs_y(counter,:) = coeffs_all(a + ncoeffs : a + 2*ncoeffs - 1);
        coeffs_z(counter,:) = coeffs_all(a + 2*ncoeffs : a + 3*ncoeffs - 1);

        counter = counter + 1;
    end
end

end
