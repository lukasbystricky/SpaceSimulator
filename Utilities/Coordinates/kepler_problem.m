function [r1, v1] = kepler_problem(r0, v0, dt)

MU = 1;
r0_sc = norm(r0);

alpha = (2*(MU/r0_sc) - norm(v0)^2)/MU;

t = @(x) (x^3*Sz(alpha*x^2) + x^2*r0'*v0*Cz(alpha*x^2)/sqrt(MU) + ...
            r0_sc*x*(1 - alpha*x^2*Sz(alpha*x^2)))/sqrt(MU);
t_dx = @(x) (x^2*Cz(alpha*x^2) + (r0'*v0)*x*(1 - alpha*x^2*Sz(alpha*x^2))/sqrt(MU) ...
            + r0_sc*(1 - alpha*x^2*Cz(alpha*x^2)))/sqrt(MU);

%%newton iteration to solve for x
TOL = 1e-8;
MAX_ITER = 50;
% x_old = alpha*dt*sqrt(mu);
x_old = 0;
res = 2*TOL;
iter = 0;
while (res > TOL && iter < MAX_ITER)
    x_new = x_old + (dt - t(x_old))/t_dx(x_old);
    res = abs(t(x_new) - t(x_old));
    
    x_old = x_new;
    iter = iter + 1;
end

if (iter == MAX_ITER)
    disp('Kepler problem failed to converge');
end

x = x_new;
z = alpha*x^2;

f = 1 - x^2*Cz(z)/r0_sc;
g = dt - x^3*Sz(z);

r1 = f*r0 + g*v0;

f_dot = sqrt(MU)/(r0_sc*norm(r1))*x*(z*Sz(z) - 1);
g_dot = 1 - x^2*Cz(z)/norm(r1);

v1 = f_dot*r0 + g_dot*v0;
    
end

function C = Cz(z)

TOL = 1e-3;

if (z > TOL)
    C = (1 - cos(sqrt(z)))/z;
else if (z < -TOL)
    C = (1 - cosh(sqrt(-z)))/z;
else
    C = 0;
    for k = 0:5
        C = C + (-z)^k/factorial(2*k + 2);
    end
    end
end
end

function S = Sz(z)
    
TOL = 1e-3;

if (z > TOL)
    S = (sqrt(z) - sin(sqrt(z)))/sqrt(z^3);
else if (z < -TOL)
    S = (sinh(sqrt(-z)) - sqrt(-z))/sqrt((-z)^3);
else
    S = 0;
    for k = 0:5
        S = S + (-z)^k/factorial(2*k + 3);
    end
    end
end
end