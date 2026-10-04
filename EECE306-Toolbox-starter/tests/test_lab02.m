% sphercial point roundtrip
r = randn(1000,3);
s = em.coord.c2sph(r);
r2 = em.coord.sph2c(s);
em.test.assertClose(max(em.vec.mag(r2-r)), 0, 1e-10, 'sph roundtrip');
% cylindrical point roundtrip
r = randn(1000,3);
c = em.coord.c2cyl(r);
r2 = em.coord.cyl2c(c);
em.test.assertClose(max(em.vec.mag(r2 - r)), 0, 1e-10,'cyl roundtrip');
% vector component roundtrip
A = randn(1000,3);
r = randn(1000,3);

As = em.coord.vecC2Sph(A, r);
A2 = em.coord.vecSph2C(As, r);

em.test.assertClose(max(em.vec.mag(A2 - A)), 0, 1e-10,'sph vector roundtrip');

A = randn(1000,3);
r = randn(1000,3);

Ac = em.coord.vecC2Cyl(A, r);
A2 = em.coord.vecCyl2C(Ac, r);

em.test.assertClose(max(em.vec.mag(A2 - A)), 0, 1e-10, 'cyl vector roundtrip');