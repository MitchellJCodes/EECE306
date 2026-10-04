function p = vecC2Cyl(A, r)

Ax = A(:,1);
Ay = A(:,2);
Az = A(:,3);

% r describes the location of vector A.
% Convert r to cylindrical coordinates to find phi.

cyl = em.coord.c2cyl(r);

phi = cyl(:,2);

A_rho = Ax .* cos(phi) + Ay .* sin(phi);

A_phi = -Ax .* sin(phi) + Ay .* cos(phi);

p = [A_rho, A_phi, Az];

end