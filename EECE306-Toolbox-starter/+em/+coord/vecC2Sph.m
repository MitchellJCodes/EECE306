function p = vecC2Sph(A, r)

Ax = A(:,1);
Ay = A(:,2);
Az = A(:,3);

sph = em.coord.c2sph(r);

theta = sph(:,2);
phi = sph(:,3);

A_r = Ax .* sin(theta) .* cos(phi) + Ay .* sin(theta) .* sin(phi) + Az .* cos(theta);

A_theta = Ax .* cos(theta) .* cos(phi) + Ay .* cos(theta) .* sin(phi) - Az .* sin(theta);

A_phi = -Ax .* sin(phi) + Ay .* cos(phi);

p = [A_r, A_theta, A_phi];

end