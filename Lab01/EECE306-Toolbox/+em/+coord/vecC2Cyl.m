function Ac = vecC2Cyl(A,r)
%VECC2CYL Summary of this function goes here
%   Detailed explanation goes here
Ax = A(1);
Ay = A(2);
z = A(3);
r_cyl = em.coord.c2cyl(r);
phi = r_cyl(2);
Arho = Ax .* cos(phi) + Ay .* sin(phi) + z .* 0;
Aphi = Ax .* (-sin(phi)) + Ay .* cos(phi) + z .* 0;
Az = z;
Ac = [Arho, Aphi, Az];
end

