%% EECE 306 Lab 5 notebook. Flux, Gauss's Law, and the Divergence
% *Team 3* Carissa McWilliams, Amelia Harris, and Mitch Nazareth
%
% Fill in every TODO, run |publish('Lab05_notebook.m')| from the toolbox
% root, print the HTML to PDF, three to five pages.

%% Setup
clear; close all;
Q = 3e-9;

%% Gauss's law on the sphere
s = em.src.pointCharge(Q, [0 0 0]);
FD = @(r) em.field.D(s, r);
sphR = @(th,ph) 0.5.*[sin(th).*cos(ph) sin(th).*sin(ph) cos(th)];
Phi = em.field.flux(FD, sphR, [0 pi], [0 2*pi], 60, 120);
fprintf('flux through sphere   = %.6e C\n', Phi);
fprintf('enclosed charge       = %.6e C\n', Q);
fprintf('relative error        = %.3e\n', abs(Phi - Q)/Q);

%% The same charge outside the surface
sout = em.src.pointCharge(Q, [5 0 0]);
Phi0 = em.field.flux(@(r) em.field.D(sout, r), sphR, [0 pi], [0 2*pi], 60, 120);
fprintf('flux with charge outside = %.3e C  (consistent with zero)\n', Phi0);
% Amelia Harris, Carissa McWilliams, and Mitch Nazareth justify the absolute tolerance you would use to call this zero.
% Compare it against the flux each hemisphere carries separately.
% We choose to go with an absoulte tolerance of 1e-12 C to consider the
% flux zero. The calculated flux of 3.379e-16 C is serval orders of
% magnitude below this tolerance and is also negligble compared with the
% flux carried by each hemisphere. Therefore, the net flux is consistent
% with zero.
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
hemi = @(th,ph) 0.5.*[sin(th).*cos(ph) sin(th).*sin(ph) cos(th)];
capD = @(u,v) [u.*cos(v) -u.*sin(v) zeros(size(u))];
PhiHemi = em.field.flux(FD3, hemi, [0 pi/2], [0 2*pi], 60, 120) ...
        + em.field.flux(FD3, capD, [0 0.5], [0 2*pi], 60, 120);
fprintf('sphere              %.6e C\n', PhiSph);
fprintf('cube                %.6e C\n', cube);
fprintf('hemisphere and cap  %.6e C\n', PhiHemi);
% The cube was the hardest surface to set up correctly because each of its
% six faces needed to have the correct outward orientation. I checked each
% face by using the order of the two surface parameters so that dr/du cross
% dr/dv pointed away from the center of the cube. If a face had the wrong
% orientation, its flux contribution would have the wrong sign.
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
    ball = @(u,v,w)[w.*sin(u).*cos(v) w.*sin(u).*sin(v) w.*cos(u)];
    sb = em.src.volCharge(rhov, ball, [0 pi], [0 2*pi], [0 0.5], Nu, 2*Nu, Nu);
    dv = em.op.div(@(r) em.field.D(sb, r), [0.1 0.05 0.1], hdiv);
    fprintf('Nu = %2d   div D inside = %.4e   rho_v = %.4e   rel err = %.2e\n', ...
        Nu, dv, rhov, abs(dv - rhov)/rhov);
end
dfree = em.op.div(FD, [2 1 1], hdiv);
fprintf('div D in charge free space = %.3e  (consistent with zero)\n', dfree);
% Carissa McWilliams, Mitch Nazareth, and Amelia Harris one or two sentences. Why must h span many source elements here,
% and what would div D return with h far smaller than the element spacing.
% The step h must span many source elements so that the finite difference
% sees the uniform charge distribution as a continuous volume rather than
% as individual discrete source elements. If h were much smaller than the
% element spacing, div D would be dominated by the discretization and could
% give an inaccurate value instead of approaching rho_v.

%% Orientation is a convention you must own
sphFlip = @(ph,th) 0.5.*[sin(th).*cos(ph) sin(th).*sin(ph) cos(th)];
PhiFlip = em.field.flux(FD, sphFlip, [0 2*pi], [0 pi], 120, 60);
fprintf('flux with swapped parameter order = %.6e C  (sign flipped)\n', PhiFlip);

%% Interpretation
% Amelia Harris, Carissa McWilliams, and Mitch Nazareth three to six sentences. Gauss's law was verified without knowing
% any closed form for the flux integral. Explain why this class of test
% remains available on a problem with no known answer.
% Gauss's law gives a physical relationship between the total flux through
% a closed surface and the charge enclosed by that surface. Therefore, I
% can test my numerical result by comparing the calculated flux with the
% known enclosed charge without needing a closed-form solution for the
% flux integral. This type of test is useful for more complicated
% geometries because the conservation law must still be satisfied. Agreement
% between different closed surfaces also gives confidence that the numerical
% implementation is correct.

%% Problems encountered
% TODO honest account, or NONE.
% The main problem encountered was getting the relative error within the
% tolerance required by the lab. Using midpoint quadrature with the
% original number of points produced results close to the expected values,
% but the relative error was still outside the required bounds. Increasing
% the number of quadrature points showed that the error decreased as the
% resolution increased. Using Gauss quadrature in the test cases provided
% enough numerical accuracy to satisfy the required relative error.

%% Full test suite
runTests
