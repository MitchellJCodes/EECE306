function test_02 = test_lab02()
%TEST_LAB02 Summary of this function goes here
%   Detailed explanation goes here
clear; close all;
rng(306);
r = randn(1000,3);
esph = max(em.vec.mag(em.coord.sph2c(em.coord.c2sph(r)) - r));
ecyl = max(em.vec.mag(em.coord.cyl2c(em.coord.c2cyl(r)) - r));
fprintf('spherical   roundtrip max error = %.3e\n', esph);
fprintf('cylindrical roundtrip max error = %.3e\n', ecyl);
disp('c2sph([0 0 1])  ='); disp(em.coord.c2sph([0 0 1]))
disp('c2sph([0 0 -1]) ='); disp(em.coord.c2sph([0 0 -1]))
q = [ 1 1 0; -1 1 0; -1 -1 0; 1 -1 0 ];
p = em.coord.c2cyl(q);
fprintf('phi in the four quadrants = %.4f %.4f %.4f %.4f rad\n', p(:,2));
fprintf('all inside [0, 2*pi) = %d\n', all(p(:,2) >= 0 & p(:,2) < 2*pi));
A  = randn(1000,3);
As = em.coord.vecC2Sph(A, r);
ev = max(em.vec.mag(em.coord.vecSph2C(As, r) - A));
em_ = max(abs(em.vec.mag(As) - em.vec.mag(A)));
fprintf('vector roundtrip (spherical pair) max error = %.3e\n', ev);
fprintf('magnitude change under conversion max       = %.3e\n', em_);
Ac = em.coord.vecC2Cyl(A, r);
ev2 = max(em.vec.mag(em.coord.vecCyl2C(Ac, r) - A));
fprintf('vector roundtrip (cylindrical pair) max error = %.3e\n', ev2);
A0 = [2 3 0];
disp('at (0, 5, 0), (Ar, Atheta, Aphi) ='); disp(em.coord.vecC2Sph(A0, [0 5 0]))
disp('at (4, 0, 0), (Ar, Atheta, Aphi) ='); disp(em.coord.vecC2Sph(A0, [4 0 0]))
end

