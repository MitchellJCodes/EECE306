function p = vecCyl2C(Ac, r)

Ac_rho = Ac(:,1);
Ac_phi = Ac(:,2);
Ac_z = Ac(:,3);

% r describes the location of vector A.
% Convert r to cylindrical coordinates to find phi.

cyl = em.coord.c2cyl(r);

phi = cyl(:,2);

Ax = Ac_rho .* cos(phi) - Ac_phi .* sin(phi);

Ay = Ac_rho .* sin(phi) + Ac_phi .* cos(phi);

p = [Ax, Ay, Ac_z];

end