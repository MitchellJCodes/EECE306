function A = vecCyl2C(Ac,r)
%VECCYL2C Summary of this function goes here
%   A = EM.COORD.VECCYL2C(Ac, r) returns the corresponding Nx3 Cartesian 
%   vector components of a set of Nx3 cylindrical vector component along an
%   Nx3 vector in the Cartesian plane.
%
%   Example:
%   em.coord.vecCyl2C([5, 0, 5], [4, -3, 2]) % returns [4, -3, 5]

Arho = Ac(:, 1);
Aphi = Ac(:, 2);
Az = Ac(:, 3);
r_cyl = em.coord.c2cyl(r);
phi = r_cyl(:, 2);
Ax = Arho .* cos(phi) + Aphi .* (-sin(phi)) + Az .* 0;
Ay = Arho .* sin(phi) + Aphi .* cos(phi) + Az .* 0;
A = [Ax, Ay, Az];
end
