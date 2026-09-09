%% EECE 306 Lab 3 notebook. Point Charges, Superposition, and the Field Engine
% *Team 3 Carissa McWilliams, Amelia Harris, and Mitch Nazareth* TODO replace with your team number and member names.
%
% Fill in every TODO, run |publish('Lab03_notebook.m')| from the toolbox
% root, print the HTML to PDF, three to five pages.

%% Setup
clear; close all;
rng(306);
e0 = em.const.eps0();

%% Closed form comparison at ten radii
Q = 2e-9;
s = em.src.pointCharge(Q, [0 0 0]);
radii = (1:10)';
Enum  = em.field.E(s, [radii zeros(10,2)]);
Eref  = Q ./ (4*pi*e0*radii.^2);
relerr = abs(Enum(:,1) - Eref) ./ Eref;
disp('   r (m)    relative error')
disp([radii relerr])

%% Symmetry on the bisector plane
% Two equal charges on the x axis. On the plane x = 0 the x component
% must vanish.
q1 = em.src.pointCharge(1e-9, [-0.1 0 0]);
q2 = em.src.pointCharge(1e-9, [ 0.1 0 0]);
s2 = em.src.merge(q1, q2);
pts = [zeros(200,1) randn(200,2)];
Ev  = em.field.E(s2, pts);
fprintf('max |Ex| on bisector      = %.3e V/m\n', max(abs(Ev(:,1))));
fprintf('max |E| on the same plane = %.3e V/m\n', max(em.vec.mag(Ev)));
% TODO justify your absolute tolerance for the first number by comparing
% it against the second. State the threshold you chose and why.
%
% The two equal charges push the electric field euqally in opposite x directions,
% which cancels out, so max |Ex| should be 0, and 10^10 should be a small enough tolerance.

%% Dipole far field
sd = em.src.merge(em.src.pointCharge(1e-9,[0 0 0.05]), em.src.pointCharge(-1e-9,[0 0 -0.05]));
E1 = em.vec.mag(em.field.E(sd, [0 0 10]));
E2 = em.vec.mag(em.field.E(sd, [0 0 20]));
fprintf('|E(r)| / |E(2r)| = %.4f   (expect near 8 for a dipole)\n', E1/E2);

%% The picture, raw and normalized
figure;
em.viz.quiver2(@(r) em.field.E(s2, r), [-0.5 0.5], [-0.5 0.5], 15, ...
    struct('normalize', false, 'title', 'Two equal charges, raw magnitudes'));
figure;
em.viz.quiver2(@(r) em.field.E(s2, r), [-0.5 0.5], [-0.5 0.5], 15, ...
    struct('normalize', true, 'title', 'Two equal charges, direction only'));
% TODO one or two sentences on why the raw plot is unreadable and what
% information the normalized plot gives up in exchange.
%
% Normalizing the plot makes it easier to see the direction of the field grid,
% as the arrows would otherwise have varying lengths. This however loses the intesity at ecah point.

%% Timing at N = 5000
pts5k = randn(5000,3) + 5;
tic; em.field.E(s2, pts5k); t = toc;
fprintf('em.field.E on 5000 points took %.3f s  (requirement, under 2 s)\n', t);

%% Interpretation
% Amelia Harris, Mitch Nazareth, and Carissa McWilliams three to six sentences. Three tests were run above, closed form,
% symmetry, and a scaling law. Say what each one can catch that the other
% two cannot.

%% Problems encountered
% Mitch Nazareth, Amelia Harris, and Carissa McWilliams honest account, or NONE.
%
% We had some mistakes from Lab 1 and Lab 2 to correct before we could begin on this document.
% Particularly, we did not specify some vectors to have Nx3 inputs. Leading to tests failing
% Interestingly enough, the tests still ran (and reported passing) when run using Octave,
% althogh we know that they were specified wrong.
%
% Finally, version controlling, with seperate branches for each individual's work on github
% massively enhanced our ability to collaborate! :)

%% Full test suite
runTests
