function test_05 = test_lab05()
%TEST_LAB05 Lab 5 checks.
%   Tests electric flux density, flux integration, and divergence.
%

clear; close all;

%% Setup

tol = 1e-10;
e0 = em.const.eps0();

%% 1. D field of a point charge

Q = 2e-9;
s = em.src.pointCharge(Q,[0 0 0]);

r = [1 0 0;
     0 2 0;
     0 0 3;
     1 2 3];

D = em.field.D(s,r);

Rmag = sqrt(sum(r.^2,2));
Dref = Q/(4*pi) .* r ./ Rmag.^3;

em.test.assertClose(D,Dref,tol, ...
    'Point charge D field');

%% 2. D = eps0 * E

E = em.field.E(s,r);

em.test.assertClose(D,e0*E,tol, ...
    'D equals eps0 times E');

%% 3. Superposition of two charges

s1 = em.src.pointCharge(1e-9,[-0.1 0 0]);
s2 = em.src.pointCharge(2e-9,[ 0.1 0 0]);
sm = em.src.merge(s1,s2);

rtest = [0.3 0.2 0.1];

Dsum = em.field.D(sm,rtest);

D1 = em.field.D(s1,rtest);
D2 = em.field.D(s2,rtest);

em.test.assertClose(Dsum,D1+D2,tol, ...
    'D field obeys superposition');

%% 4. Observation point too close to source

caught = false;

try
    em.field.D(s,[0 0 0]);
catch
    caught = true;
end

assert(caught, ...
    'D must reject observation points at source locations');

%% 5. Flux of a constant field through a unit square

Fconst = @(r) [zeros(size(r,1),2) ones(size(r,1),1)];

square = @(u,v) [u v zeros(size(u))];

Phi = em.field.flux( ...
    Fconst,square,[0 1],[0 1],20,20);

em.test.assertClose(Phi,1,1e-10, ...
    'Constant field through unit square');

%% 6. Reversing the surface orientation flips the flux

PhiFlip = em.field.flux( ...
    Fconst,square,[1 0],[0 1],20,20);

em.test.assertClose(PhiFlip,-1,1e-10, ...
    'Reversing surface orientation flips flux');

%% 7. Gauss law for a point charge

sphere = @(th,ph) ...
    0.5*[sin(th).*cos(ph) ...
         sin(th).*sin(ph) ...
         cos(th)];

PhiQ = em.field.flux( ...
    @(r) em.field.D(s,r), ...
    sphere,[0 pi],[0 2*pi],60,120);

em.test.assertClose(PhiQ,Q,2e-4, ...
    'Gauss law for point charge');

%% 8. Charge outside closed surface gives zero flux

sout = em.src.pointCharge(Q,[5 0 0]);

PhiOutside = em.field.flux( ...
    @(r) em.field.D(sout,r), ...
    sphere,[0 pi],[0 2*pi],60,120);

em.test.assertClose(PhiOutside,0,1e-10, ...
    'Flux from external charge');

%% 9. Divergence of a linear vector field

% F = [2x, 3y, 4z]
% div(F) = 2 + 3 + 4 = 9

Flinear = @(r) ...
    [2*r(:,1) 3*r(:,2) 4*r(:,3)];

r0 = [0.37 -0.21 0.44];

divF = em.op.div(Flinear,r0,1e-4);

em.test.assertClose(divF,9,1e-10, ...
    'Divergence of linear field');

%% 10. Divergence of a constant vector field

Fzero = @(r) ...
    repmat([3 -2 5],size(r,1),1);

divZero = em.op.div(Fzero,[0.2 0.3 0.4],1e-4);

em.test.assertClose(divZero,0,1e-10, ...
    'Divergence of constant field');

%% 11. Divergence of D inside a uniform charged volume

rhoV = 1e-6;
h = 0.05;

ball = @(u,v,w) ...
    [w.*sin(u).*cos(v) ...
     w.*sin(u).*sin(v) ...
     w.*cos(u)];

sb = em.src.volCharge( ...
    rhoV,ball, ...
    [0 pi],[0 2*pi],[0 0.5], ...
    40,80,40);

divD = em.op.div( ...
    @(r) em.field.D(sb,r), ...
    [0.1 0.05 0.1],h);

fprintf('div D = %.6e, rho_v = %.6e\n',divD,rhoV);

assert(abs(divD-rhoV)/rhoV < 0.1, ...
    'Divergence of D should approach the volume charge density');

disp('  test_lab05 checks complete');

end
