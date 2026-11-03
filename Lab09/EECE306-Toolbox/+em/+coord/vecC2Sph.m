function As = vecC2Sph(A,r)
%VECC2SPH Cartesian to Spherical Vector Component Conversion
%   AS = EM.COORD.VECC2SPH(A, r) returns the corresponding Nx3 spherical 
%   vector components of an Nx3 Cartesian component set along an Nx3 vector
%   in the Cartesian plane.
%
%   Example:
%   em.coord.vecC2Sph([pi/6, pi/4, 1], [0, 0, 1]) % returns [1.0000, 0.5236, 0.7854]
Ax = A(:,1);
Ay = A(:,2);
Az = A(:,3);
r_sph = em.coord.c2sph(r);
theta = r_sph(:,2);
phi = r_sph(:,3);
Aradial = Ax .* sin(theta) .* cos(phi) + Ay .* sin(theta) .* sin(phi) + Az .* cos(theta);
Atheta = Ax .* cos(theta) .* cos(phi) + Ay .* cos(theta) .* sin(phi) + Az .* (-sin(theta));
Aphi = Ax .* (-sin(phi)) + Ay .* cos(phi) + Az .* 0;
As = [Aradial, Atheta, Aphi];
end

