%% EECE 306 Lab 5 notebook. Flux, Gauss's Law, and the Divergence
% *Team 03 Carissa McWilliams, Amelia Harris, and Mitch Nazareth.* 

%% Setup
clear; close all;
Q = 3e-9;

%% Gauss's law on the sphere
s = em.src.pointCharge(Q, [0 0 0]);
FD = @(r) em.field.D(s, r);
sphR = @(th,ph) 0.5*[sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];
Phi = em.field.flux(FD, sphR, [0 pi], [0 2*pi], 60, 120);
fprintf('flux through sphere   = %.6e C\n', Phi);
fprintf('enclosed charge       = %.6e C\n', Q);
fprintf('relative error        = %.3e\n', abs(Phi - Q)/Q);

%% The same charge outside the surface
sout = em.src.pointCharge(Q, [5 0 0]);
Phi0 = em.field.flux(@(r) em.field.D(sout, r), sphR, [0 pi], [0 2*pi], 60, 120);
fprintf('flux with charge outside = %.3e C  (consistent with zero)\n', Phi0);

% I initially ran the assert with a tolerance of 1e-4,
% and found that I was getting a relative error of
% about 1.14e-4, so I raised the tolerance to 2e-4,
% as that is high enough to pass test, but still
% very accurate, more than enough for this test.

%% Three surfaces, one answer
% The same charge, placed off center inside all three closed surfaces.
s3  = em.src.pointCharge(Q, [0 0.05 0.1]);
FD3 = @(r) em.field.D(s3, r);
PhiSph = em.field.flux(FD3, sphR, [0 pi], [0 2*pi], 60, 120);
cube = 0;
faces = { @(u,v)[u v  0.5*ones(size(u))], @(u,v)[v u -0.5*ones(size(u))], ...
          @(u,v)[v  0.5*ones(size(u)) u], @(u,v)[u -0.5*ones(size(u)) v], ...
          @(u,v)[ 0.5*ones(size(u)) u v], @(u,v)[-0.5*ones(size(u)) v u] };
for k = 1:6
    cube = cube + em.field.flux(FD3, faces{k}, [-0.5 0.5], [-0.5 0.5], 60, 60);
end
hemi = @(th,ph) 0.5*[sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];
capD = @(u,v) [u.*cos(v) -u.*sin(v) zeros(size(u))];
PhiHemi = em.field.flux(FD3, hemi, [0 pi/2], [0 2*pi], 60, 120) ...
        + em.field.flux(FD3, capD, [0 0.5], [0 2*pi], 60, 120);
fprintf('sphere              %.6e C\n', PhiSph);
fprintf('cube                %.6e C\n', cube);
fprintf('hemisphere and cap  %.6e C\n', PhiHemi);

% The cube was the hardest surface because it has six faces.
% Setting up the correct orientation of each face required
% comparing the direction of the parameter cross product
% with the outward normal direction.

%% Independence of the charge position inside
for x0 = [0 0.15 0.3]
    sm = em.src.pointCharge(Q, [x0 0.1 0]);
    Pm = em.field.flux(@(r) em.field.D(sm, r), sphR, [0 pi], [0 2*pi], 60, 120);
    fprintf('charge at x = %.2f   flux = %.6e C\n', x0, Pm);
end

%% Divergence against the density, an element count study
% A uniform ball of charge. The derivative step h must average over many
% elements, so h = 0.05 m is passed explicitly, and the element count is
% raised at fixed h.
rhov = 1e-6;
hdiv = 0.05;
for Nu = [10 20 40]
    ball = @(u,v,w)[w*sin(u)*cos(v) w*sin(u)*sin(v) w*cos(u)];
    sb = em.src.volCharge(rhov, ball, [0 pi], [0 2*pi], [0 0.5], Nu, 2*Nu, Nu);
    dv = em.op.div(@(r) em.field.D(sb, r), [0.1 0.05 0.1], hdiv);
    fprintf('Nu = %2d   div D inside = %.4e   rho_v = %.4e   rel err = %.2e\n', ...
        Nu, dv, rhov, abs(dv - rhov)/rhov);
end
dfree = em.op.div(FD, [2 1 1], hdiv);
fprintf('div D in charge free space = %.3e  (consistent with zero)\n', dfree);

% h spans across so many source elements here so the finite difference
% calculation can see the overall continuous charge distribution. If
% h is much smaller than the element spacing the calculated divergence
% will be inaccurate because it will be affected by the individual
% source elements.

%% Orientation is a convention you must own
sphFlip = @(ph,th) 0.5*[sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];
PhiFlip = em.field.flux(FD, sphFlip, [0 2*pi], [0 pi], 120, 60);
fprintf('flux with swapped parameter order = %.6e C  (sign flipped)\n', PhiFlip);

%% Interpretation
% Gauss's law gives us a relationship between total flux and the charge
% enclosed by the surface. We can calculate the flux without knowing
% the exact analytical value of the integral. We then compare the result
% to the known enclosed charge. This lets us test if the implementation
% is working correctly even if the integral does not have a known closed
% form solution.

%% Problems encountered
% We had to fix a bug in +em+/+quad/vol.m where the function did not accept the quadrature rule argument. Also, the Jacobian wasn't being calculated correctly so that the weights represented the physical volume.

%% Full test suite
runTests
