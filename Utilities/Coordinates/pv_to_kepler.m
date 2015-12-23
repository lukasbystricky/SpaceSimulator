function parameters = pv_to_kepler(r, v)

h = cross(r, v);
Omega = atan2(h(1),-h(2));
i = atan2(sqrt(h(1)^2 + h(2)^2),h(3));

R1 = [cos(Omega),  sin(Omega), 0;
      -sin(Omega), cos(Omega), 0;
      0,            0,           1];

R2 = [1, 0 ,       0;
      0, cos(i),  sin(i);
      0, -sin(i), cos(i)];

  
p = R2*R1*r;
u = atan2(p(2),p(1));

a = norm(r)/(2 - norm(r)*norm(v)^2);
e = sqrt(1 - norm(h)^2/a);

cosE = (a - norm(r))/(a*e);
sinE = r'*v/(e*sqrt(a));

nu = atan2(sqrt(1 - e^2)*sinE,(cosE - e));
omega = u - nu;

parameters.a = a;
parameters.e = e;
parameters.i = i*180/pi;
parameters.omega = wrapTo2Pi(omega)*180/pi;
parameters.Omega = wrapTo2Pi(Omega)*180/pi;
parameters.nu = wrapTo2Pi(nu)*180/pi;
