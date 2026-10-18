%% EECE 306 Lab 2 notebook. Coordinate Systems, Points and Vector Components
% Team 3: Amelia Harris and Carissa McWilliams
%
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
% TODO two to four sentences. The vector did not change. Explain exactly
% what did change, and why a position argument is therefore mandatory for
% component conversion and meaningless for point conversion.
%
% The directions of the spherical components (r_hat, theta_hat, phi_hat)
% changed at each position along the vector, since they depend on the polar
% angle theta and the azimuthal angle phi. Spherical axes are local to a 
% point and rotate as the point moves through space. Thus, component 
% conversion requires a position argument to determine the orientation of
% the local axes in 3D space. Point conversion does not change the
% orientation of the axes; it only changes how the location of the point is
% described. Therefore, a position argument is unnecessary in this case.

%% Full test suite
runTests

%% Interpretation
% TODO what class of bug does the roundtrip test catch, and what class
% does it miss. Name the invariant that catches what the roundtrip
% misses, and say why it is reference free.
%
% 1. The roundtrip test catches the class of bug in which a function treats
% a cylindrical coordinate (rho, phi, z) as a Cartesian coordinate (x, y,
% z). The roundtrip test misses the class of bug in which atan2 returns
% a negative phi when x is negative.
% 2. The invariant that catches what the roundtrip misses is the magnitude
% of the vector. It is reference-free because it is unchanged by any
% component conversion and can always be calculated if the scalar
% components are known.

%% Problems encountered
% TODO honest account, or NONE.
% Amelia: Minor syntax and labeling errors when copying similar code between
% toolkit files. em.coord functions return Nx1 vector outputs, not Nx3.
% Roundtrip errors from conversions were considerably high. For test 6, phi
% only shows as one value, not four as intended. Most of these issues were
% fixed by correcting the way the function takes the individual coordinates
% from the input.
% Carissa: Not expressing point and vector conversions as empty Nx3 
% arrays and as separate or A values, which stopped the random 1000 x 3
% test the eight functions test from working properly. No single output for
% the point conversion function which lead to another error which 
% interfered with matrix multiplication.
%% Full test suite
runTests
