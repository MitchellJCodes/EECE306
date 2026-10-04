function p = vecSph2C(As, r)

As_r = As(:,1);
As_theta = As(:,2);
As_phi = As(:,3);

sph = em.coord.c2sph(r);

theta = sph(:,2);
phi = sph(:,3);

Ax = As_r .* sin(theta) .* cos(phi) + As_theta .* cos(theta) .* cos(phi) - As_phi .* sin(phi);

Ay = As_r .* sin(theta) .* sin(phi) + As_theta .* cos(theta) .* sin(phi) + As_phi .* cos(phi);

Az = As_r .* cos(theta) - As_theta .* sin(theta);

p = [Ax, Ay, Az];

end