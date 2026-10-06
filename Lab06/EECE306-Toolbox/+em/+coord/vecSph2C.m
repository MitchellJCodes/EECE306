function A = vecSph2C(As, r)
%VECSPH2C Spherical to Cartesian Vector Component Conversion
%   AC = EM.COORD.VECSPH2C(A, r) returns the corresponding Nx3 Cartesian
%   vector components  of a set of Nx3 spherical vector component along an
%   Nx3 vector in the Cartesian plane.
%
%   Example:
%   em.coord.vecSph2C([1, pi/6, pi/4], [0, 0, 1]) % returns [0.5236, 0.7854, 1.0000]
Aradial = As(:, 1);
Atheta = As(:, 2);
Aphi = As(:, 3);
r_sph = em.coord.c2sph(r);
theta = r_sph(:, 2);
phi = r_sph(:, 3);
Ax = Aradial .* sin(theta) .* cos(phi) + Atheta .* cos(theta) .* cos(phi) + Aphi .* (-sin(phi));
Ay = Aradial .* sin(theta) .* sin(phi) + Atheta .* cos(theta) .* sin(phi) + Aphi .* cos(phi);
Az = Aradial .* cos(theta) + Atheta .* (-sin(theta)) + Aphi .* 0;
A = [Ax, Ay, Az];
end

