r = randn(1000,3);
s = em.coord.c2sph(r);
r2 = em.coord.sph2c(s);
em.test.assertClose(max(em.vec.mag(r2-r)), 0, 1e-10, 'sph roundtrip');