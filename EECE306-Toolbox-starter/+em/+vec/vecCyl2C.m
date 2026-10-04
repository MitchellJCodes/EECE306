function p = vecCyl2C(Ac, r)
Ac_rho = Ac(1)
Ac_phi = Ac(2);
Ac(z) = Ac(3);

% r describes the location of vector A.
% Convert r to cylindrical coordinates to find phi.
[rho, phi, z] = c2cyl(r);
 Ax = Ac_rho .* cos(phi) - Ac_phi .* sin(phi);
 Ay = rho * sin(phi);

 p = [Ax, Ay, z]
end