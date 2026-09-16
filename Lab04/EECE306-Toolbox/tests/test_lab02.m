function test_02 = test_lab02()
%TEST_LAB02 Tests coordinate conversion functions.
%   Generates random Nx3 Cartesian points and
%   then verifies conversions to spherical
%   and cylindrical coordinate systems.
%   Checks spherical pole values and cylinder
%   angles in all 4 quadrants. Converts vector
%   components and checks for correct magnitudes.
%   Re-written to work with assertions.
clear; close all;
rng(306);

tol = 1e-10;

r = randn(1000,3);

% Cartesian -> spherical -> Cartesian
esph = max(em.vec.mag(em.coord.sph2c(em.coord.c2sph(r)) - r));
fprintf('spherical   roundtrip max error = %.3e\n', esph);
em.tests.assertClose(esph, 0, tol, 'Spherical roundtrip error');

% Cartesian -> cylindrical -> Cartesian
ecyl = max(em.vec.mag(em.coord.cyl2c(em.coord.c2cyl(r)) - r));
fprintf('cylindrical roundtrip max error = %.3e\n', ecyl);
em.tests.assertClose(ecyl, 0, tol, 'Cylindrical roundtrip error');

disp('c2sph([0 0 1])  ='); 
disp(em.coord.c2sph([0 0 1]))

disp('c2sph([0 0 -1]) ='); 
disp(em.coord.c2sph([0 0 -1]))

q = [ 1 1 0; -1 1 0; -1 -1 0; 1 -1 0 ];
p = em.coord.c2cyl(q);

fprintf('phi in the four quadrants = %.4f %.4f %.4f %.4f rad\n', p(:,2));

expected_phi = [pi/4; 3*pi/4; 5*pi/4; 7*pi/4];
em.tests.assertClose(p(:,2), expected_phi, tol, ...
    'Phi values in the four quadrants');

inside_range = all(p(:,2) >= 0 & p(:,2) < 2*pi);
fprintf('all inside [0, 2*pi) = %d\n', inside_range);
em.tests.assertClose(inside_range, 1, tol, ...
    'Phi range');

A  = randn(1000,3);

% Cartesian -> spherical vector components -> Cartesian
As = em.coord.vecC2Sph(A, r);
ev = max(em.vec.mag(em.coord.vecSph2C(As, r) - A));
fprintf('vector roundtrip (spherical pair) max error = %.3e\n', ev);
em.tests.assertClose(ev, 0, tol, ...
    'Spherical vector roundtrip');

% Check that vector magnitude is preserved
em_ = max(abs(em.vec.mag(As) - em.vec.mag(A)));
fprintf('magnitude change under conversion max       = %.3e\n', em_);
em.tests.assertClose(em_, 0, tol, ...
    'Spherical vector magnitude');

% Cartesian -> cylindrical vector components -> Cartesian
Ac = em.coord.vecC2Cyl(A, r);
ev2 = max(em.vec.mag(em.coord.vecCyl2C(Ac, r) - A));
fprintf('vector roundtrip (cylindrical pair) max error = %.3e\n', ev2);
em.tests.assertClose(ev2, 0, tol, ...
    'Cylindrical vector roundtrip');

A0 = [2 3 0];

disp('at (0, 5, 0), (Ar, Atheta, Aphi) ='); 
disp(em.coord.vecC2Sph(A0, [0 5 0]))

disp('at (4, 0, 0), (Ar, Atheta, Aphi) ='); 
disp(em.coord.vecC2Sph(A0, [4 0 0]))

end
