function A = vecCyl2C(Ac,r)
%VECC2CYL Summary of this function goes here
%   Detailed explanation goes here
Arho = Ac(1);
Aphi = Ac(2);
z = Ac(3);
r_cyl = em.coord.c2cyl(r);
phi = r_cyl(2);
Ax = Arho .* cos(phi) + Aphi .* (-sin(phi)) + z .* 0;
Ay = Arho .* sin(phi) + Aphi .* cos(phi) + z .* 0;
Az = z;
A = [Ax, Ay, Az];
end
