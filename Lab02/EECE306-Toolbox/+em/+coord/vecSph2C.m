function A = vecSph2C(As, r)
%VECSPH2C Summary of this function goes here
%   Detailed explanation goes here
Aradial = As(1);
Atheta = As(2);
Aphi = As(3);
r_sph = em.coord.c2sph(r);
theta = r_sph(2);
phi = r_sph(3);
Ax = Aradial .* sin(theta) .* cos(phi) + Atheta .* cos(theta) .* cos(phi) + Aphi .* (-sin(phi));
Ay = Aradial .* sin(theta) .* sin(phi) + Atheta .* cos(theta) .* sin(phi) + Aphi .* cos(phi);
Az = Aradial .* cos(theta) + Atheta .* (-sin(theta)) + Aphi .* 0;
A = [Ax, Ay, Az];
end

