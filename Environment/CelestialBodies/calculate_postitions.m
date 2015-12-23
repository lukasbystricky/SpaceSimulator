function pos = calculate_postitions(ephemeris, t, planet)

switch planet
    case 'mercury'
        time_subintervals = ephemeris.mercury.time_subintervals;
        coeffs_x = ephemeris.mercury.coeffs_x;
        coeffs_y = ephemeris.mercury.coeffs_y;
        coeffs_z = ephemeris.mercury.coeffs_z;
        
    case 'venus'
        time_subintervals = ephemeris.venus.time_subintervals;
        coeffs_x = ephemeris.venus.coeffs_x;
        coeffs_y = ephemeris.venus.coeffs_y;
        coeffs_z = ephemeris.venus.coeffs_z;
        
    case 'earth'
        time_subintervals = ephemeris.earth.time_subintervals;
        coeffs_x = ephemeris.earth.coeffs_x;
        coeffs_y = ephemeris.earth.coeffs_y;
        coeffs_z = ephemeris.earth.coeffs_z;
        
    case 'mars'
        time_subintervals = ephemeris.mars.time_subintervals;
        coeffs_x = ephemeris.mars.coeffs_x;
        coeffs_y = ephemeris.mars.coeffs_y;
        coeffs_z = ephemeris.mars.coeffs_z;
        
    case 'jupiter'
        time_subintervals = ephemeris.jupiter.time_subintervals;
        coeffs_x = ephemeris.jupiter.coeffs_x;
        coeffs_y = ephemeris.jupiter.coeffs_y;
        coeffs_z = ephemeris.jupiter.coeffs_z;
        
    case 'saturn'
        time_subintervals = ephemeris.saturn.time_subintervals;
        coeffs_x = ephemeris.saturn.coeffs_x;
        coeffs_y = ephemeris.saturn.coeffs_y;
        coeffs_z = ephemeris.saturn.coeffs_z;
        
    case 'uranus'
        time_subintervals = ephemeris.uranus.time_subintervals;
        coeffs_x = ephemeris.uranus.coeffs_x;
        coeffs_y = ephemeris.uranus.coeffs_y;
        coeffs_z = ephemeris.uranus.coeffs_z;
        
    case 'neptune'
        time_subintervals = ephemeris.neptune.time_subintervals;
        coeffs_x = ephemeris.neptune.coeffs_x;
        coeffs_y = ephemeris.neptune.coeffs_y;
        coeffs_z = ephemeris.neptune.coeffs_z;
        
    case 'pluto'
        time_subintervals = ephemeris.pluto.time_subintervals;
        coeffs_x = ephemeris.pluto.coeffs_x;
        coeffs_y = ephemeris.pluto.coeffs_y;
        coeffs_z = ephemeris.pluto.coeffs_z;
        
    case 'moon'
        time_subintervals = ephemeris.moon.time_subintervals;
        coeffs_x = ephemeris.moon.coeffs_x;
        coeffs_y = ephemeris.moon.coeffs_y;
        coeffs_z = ephemeris.moon.coeffs_z;
        
    case 'sun'
        time_subintervals = ephemeris.sun.time_subintervals;
        coeffs_x = ephemeris.sun.coeffs_x;
        coeffs_y = ephemeris.sun.coeffs_y;
        coeffs_z = ephemeris.sun.coeffs_z;
end

pos = zeros(length(t), 3);

[n_subintervals, ~] = size(coeffs_x);

for i = 1:length(t)
   
   if (t(i) < time_subintervals(1) || t(i) > time_subintervals(end))
       disp('time out of range');
       break;
   else
       %find time subinterval
       for j = 1:n_subintervals
           if (t(i) >= time_subintervals(j,1) && t(i) <= time_subintervals(j,2))
               subinterval = j;
           end
       end
   end
   
   a = time_subintervals(subinterval,1);
   b = time_subintervals(subinterval,2);
   
   scaled_time = (t(i) - (a + b)/2)*(2/(b - a));
   pos(i,1) = chebyshev_approx(scaled_time, coeffs_x(subinterval,:));
   pos(i,2) = chebyshev_approx(scaled_time, coeffs_y(subinterval,:));
   pos(i,3) = chebyshev_approx(scaled_time, coeffs_z(subinterval,:));
   
end


end

function p = chebyshev_approx(t, coeffs)

p = 0;

for i = 0:length(coeffs)-1
   p = p + coeffs(i+1)*chebyshev_polynomial(t, i); 
end
end
