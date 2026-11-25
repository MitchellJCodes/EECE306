%% EECE 306 Lab 2 notebook. Coordinate Systems, Points and Vector Components
% *Team 3.* Carissa McWilliams and Ameila Harris

% Fill in every TODO, run |publish('Lab02_notebook.m')| from the toolbox
% root, print the HTML to PDF, three to five pages.

%% Setup
clear; close all;
rng(306);

%% Point roundtrips, 1000 random points
r = randn(1000,3);
esph = max(em.vec.mag(em.coord.sph2c(em.coord.c2sph(r)) - r));
ecyl = max(em.vec.mag(em.coord.cyl2c(em.coord.c2cyl(r)) - r));
fprintf('spherical   roundtrip max error = %.3e\n', esph);
fprintf('cylindrical roundtrip max error = %.3e\n', ecyl);

%% Poles and the phi convention
disp('c2sph([0 0 1])  ='); disp(em.coord.c2sph([0 0 1]))
disp('c2sph([0 0 -1]) ='); disp(em.coord.c2sph([0 0 -1]))
q = [ 1 1 0; -1 1 0; -1 -1 0; 1 -1 0 ];
p = em.coord.c2cyl(q);
fprintf('phi in the four quadrants = %.4f %.4f %.4f %.4f rad\n', p(:,2));
fprintf('all inside [0, 2*pi) = %d\n', all(p(:,2) >= 0 & p(:,2) < 2*pi));

%% Vector component roundtrips and the invariant
A  = randn(1000,3);
As = em.coord.vecC2Sph(A, r);
ev = max(em.vec.mag(em.coord.vecSph2C(As, r) - A));
em_ = max(abs(em.vec.mag(As) - em.vec.mag(A)));
fprintf('vector roundtrip (spherical pair) max error = %.3e\n', ev);
fprintf('magnitude change under conversion max       = %.3e\n', em_);
Ac = em.coord.vecC2Cyl(A, r);
ev2 = max(em.vec.mag(em.coord.vecCyl2C(Ac, r) - A));
fprintf('vector roundtrip (cylindrical pair) max error = %.3e\n', ev2);

%% The two position demonstration
% The same Cartesian vector A = 2x + 3y, expressed in spherical
% components at two different positions.
A0 = [2 3 0];
disp('at (0, 5, 0), (Ar, Atheta, Aphi) ='); disp(em.coord.vecC2Sph(A0, [0 5 0]))
disp('at (4, 0, 0), (Ar, Atheta, Aphi) ='); disp(em.coord.vecC2Sph(A0, [4 0 0]))

% The components change even with the same A0 becuase the r vectors are at
% different locations. When expressing A0 in spherical coordinates the
% source vector changes position of the converted coordinates. When you
% convert back to cartesian plane the two positions convert back to A0.



%% Interpretation
% Carissa McWilliams and Amelia Harris what class of bug does the roundtrip test catch, and what class
% does it miss. Name the invariant that catches what the roundtrip
% misses, and say why it is reference free.

%% Problems encountered
% Carissa McWilliams and Amelia Harris honest account, or NONE.
% Some common errors that appeared in tests where not expressing my point
% and vector conversions as empty N x3 arrays and as separate or A
% values. This affect my random 1000 x 3 test which restricted the program
% from testing of the eight functions works properly. I also didn't have a
% single output for my point conversion function which lead to another error
% which interfered with my matrix multiplication.
%% Full test suite
runTests
