q1 = em.src.pointCharge(1e-9, [-0.1 0 0]);
q2 = em.src.pointCharge(1e-9, [ 0.1 0 0]);
s = em.src.merge(q1, q2);
p = [zeros(200,1), randn(200,1), randn(200,1)]; % the x = 0 plane
Ev = em.field.E(s, p);
tol = 1e-6 * max(em.vec.mag(Ev)); % tiny compared with the field itself
em.test.assertClose(max(abs(Ev(:,1))), 0, tol, 'bisector symmetry');