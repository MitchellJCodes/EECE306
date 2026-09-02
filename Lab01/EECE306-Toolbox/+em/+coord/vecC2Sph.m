function As = vecC2Sph(A,r)
%VECC2SPH Summary of this function goes here
%   Detailed explanation goes here
Ax = A(1);
Ay = A(2);
Az = A(3);
r_sph = em.coord.c2sph(r);
theta = r_sph(2);
phi = r_sph(3);
Aradial = Ax .* sin(theta) .* cos(phi) + Ay .* sin(theta) .* sin(phi) + Az .* cos(theta);
Atheta = Ax .* cos(theta) .* cos(phi) + Ay .* cos(theta) .* sin(phi) + Az .* (-sin(theta));
Aphi = Ax .* (-sin(phi)) + Ay .* cos(phi) + Az .* 0;
As = [Aradial, Atheta, Aphi];
end

