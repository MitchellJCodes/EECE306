function Ac = vecC2Cyl(A,r)
%VECC2CYL Cartesian to Cylindrical Vector Component Conversion
%   AC = EM.COORD.VECC2CYL(A, r) returns the corresponding Nx3 cylindrical 
%   vector components of an Nx3 Cartesian component set along an Nx3 vector
%   in the Cartesian plane.
%
%   Example:
%   em.coord.vecC2Cyl([4, -3, 5], [4, -3, 2]) % returns [5, 0, 5]

Ax = A(:,1);
Ay = A(:,2);
Az = A(:,3);
r_cyl = em.coord.c2cyl(r);
phi = r_cyl(:,2);
Arho = Ax .* cos(phi) + Ay .* sin(phi) + Az .* 0;
Aphi = Ax .* (-sin(phi)) + Ay .* cos(phi) + Az .* 0;
Ac = [Arho, Aphi, Az];
end

