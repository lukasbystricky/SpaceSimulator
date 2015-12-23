function nu = true_anomaly_from_mean(M, ecc)

MAX_ITER = 10;
TOL = 1e-8;
E = newton_raphson(@(E) E - ecc*sin(E) - M, @(E) 1 - ecc*cos(E), M, MAX_ITER, TOL);

nu = wrapTo2Pi(2*atan2(sqrt(1 + ecc)*sin(E/2), sqrt(1 - ecc)*cos(E/2)));

