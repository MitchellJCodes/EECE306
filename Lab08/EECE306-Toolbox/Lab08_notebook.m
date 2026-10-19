%% EECE 306 Lab 8 notebook. Dielectrics, Boundary Conditions, and Capacitance
% *Team 03 — Carissa McWilliams, Amelia Harris, and Mitch Nazareth.*
% Milestone week. The graded deliverable is the formal report, which
% replaces the weekly notebook. This script generates every
% number and figure of the validation ladder so the report can cite
% them. Run |publish('Lab08_notebook.m')| and lift the results.

%% Setup
clear; close all;
e0 = em.const.eps0();

%% Rung 1. Parallel plate against the formula
% Plates span the full width, so the field between them is exactly one
% dimensional and C per unit depth must equal eps0 times w over d.
n = 81;
G = em.solve.grid2d([0 1], [0 1], n, n);
iy1 = round(0.25*(n-1)) + 1;  iy2 = round(0.75*(n-1)) + 1;
mask = false(n,n); mask(iy1,:) = true; mask(iy2,:) = true;
Vfix = zeros(n,n); Vfix(iy2,:) = 1;
C = em.derive.capacitance(G, 1, mask, Vfix, struct('method','direct'));
Cref = e0 * 1 / 0.5;
fprintf('parallel plate  C = %.5e   formula %.5e   rel err %.2e\n', C, Cref, abs(C-Cref)/Cref);

%% Rung 2. Two dielectric layers in series
epsr = ones(n,n); epsr(G.Y < 0.5) = 4;
C2 = em.derive.capacitance(G, epsr, mask, Vfix, struct('method','direct'));
C2ref = e0 / (0.25/4 + 0.25/1);
fprintf('two layers      C = %.5e   series formula %.5e   rel err %.2e\n', C2, C2ref, abs(C2-C2ref)/C2ref);

%% Rung 3. Coaxial line against the log formula
n = 121;
Gc = em.solve.grid2d([-1 1], [-1 1], n, n);
R = sqrt(Gc.X.^2 + Gc.Y.^2);
a = 0.3; b = 0.9;
maskC = R <= a | R >= b;
VfixC = zeros(n,n); VfixC(R <= a) = 1;
Cc = em.derive.capacitance(Gc, 1, maskC, VfixC, struct('method','direct'));
CcRef = 2*pi*e0 / log(b/a);
fprintf('coax            C = %.5e   formula %.5e   rel err %.2e\n', Cc, CcRef, abs(Cc-CcRef)/CcRef);
% TODO one sentence. The boundary is staircased on a Cartesian grid,
% state which direction the error goes as the grid is refined.

%% The two blunder catchers
VfixC100 = zeros(n,n); VfixC100(R <= a) = 100;
Cc100 = em.derive.capacitance(Gc, 1, maskC, VfixC100, struct('method','direct'));
fprintf('C at V0 = 1     %.5e\nC at V0 = 100   %.5e   (a missing square would differ by 1e4)\n', Cc, Cc100);
Cfill = em.derive.capacitance(Gc, 3, maskC, VfixC, struct('method','direct'));
fprintf('filled with epsr = 3, ratio C_filled over C = %.4f  (expect 3)\n', Cfill/Cc);

%% Resistance of a uniform bar
n = 81;
Gr = em.solve.grid2d([0 1], [0 1], n, n);
maskR = false(n,n); maskR(:,1) = true; maskR(:,end) = true;
VfixR = zeros(n,n); VfixR(:,1) = 1;
sigma = 100;
Rbar = em.derive.resistance(Gr, sigma, maskR, VfixR, struct('method','direct'));
Rref = 1 / (sigma * 1);
fprintf('bar resistance  R = %.5e   formula L over sigma w = %.5e\n', Rbar, Rref);

%% The payoff. A microstrip design curve
% A strip above a grounded plane with a dielectric substrate. No closed
% form exists. Every rung below this one is what authorizes the numbers.
n = 121;
Gm = em.solve.grid2d([-1 1], [0 1], n, round((n+1)/2));
ny = size(Gm.X, 1);
hsub = 0.2;
epsm = ones(ny, n); epsm(Gm.Y < hsub) = 4;
iy = find(Gm.y >= hsub, 1);
widths = [0.1 0.2 0.3 0.4 0.6];
Cw = zeros(size(widths));
for k = 1:numel(widths)
    w = widths(k);
    mk = false(ny, n);
    mk(1,:) = true; mk(end,:) = true; mk(:,1) = true; mk(:,end) = true;
    mk(iy, abs(Gm.x) <= w/2) = true;
    Vf = zeros(ny, n);
    Vf(iy, abs(Gm.x) <= w/2) = 1;
    Cw(k) = em.derive.capacitance(Gm, epsm, mk, Vf, struct('method','direct'));
end
figure; plot(widths/hsub, Cw*1e12, 'o-'); grid on
xlabel('w / h'); ylabel('C  (pF per m)')
title('Microstrip capacitance design curve, substrate epsr = 4')
disp([widths(:)/hsub Cw(:)*1e12])
% TODO two or three sentences for the report. State the trend, compare
% the widest case against the parallel plate estimate epsr eps0 w over h,
% and attribute the excess to fringing.

%% Full regression suite, Lab 1 through Lab 8
runTests
