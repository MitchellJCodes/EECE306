function test_05 = test_lab05()
%TEST_LAB05 Lab 5 Test Function
%   Runs required tests for Lab 5

%% Setup
clear; close all;
Q = 3e-9;

%% Gauss's law on the sphere
s = em.src.pointCharge(Q, [0 0 0]);
FD = @(r) em.field.D(s, r);
sphR = @(th,ph) 0.5*[sin(th).*cos(ph) sin(th).*sin(ph) cos(th)];
Phi = em.field.flux(FD, sphR, [0 pi], [0 2*pi], 60, 120, 'simpson');
fprintf('flux through sphere   = %.6e C\n', Phi);
fprintf('enclosed charge       = %.6e C\n', Q);
fprintf('relative error        = %.3e \n', abs(Phi - Q)/Q);
em.test.assertClose(Phi, 0, 1e-6, 'bisector symmetry');

%% The same charge outside the surface
sout = em.src.pointCharge(Q, [5 0 0]);
Phi0 = em.field.flux(@(r) em.field.D(sout, r), sphR, [0 pi], [0 2*pi], ...
    60, 120, 'simpson');
fprintf('flux with charge outside = %.3f C  (consistent with zero)\n', ...
    Phi0);
em.test.assertClose(Phi0, 0, 1e-6, ['Zero flux through sphere with' ...
    ' outside charge'])

%% Three surfaces, one answer
% The same charge, placed off center inside all three closed surfaces.
s3  = em.src.pointCharge(Q, [0 0.05 0.1]);
FD3 = @(r) em.field.D(s3, r);
PhiSph = em.field.flux(FD3, sphR, [0 pi], [0 2*pi], 60, 120, 'simpson');
cube = 0;
faces = { @(u,v)[u v  0.5.*ones(size(u))], @(u,v)[v u -0.5.*ones(size(u))], ...
          @(u,v)[v  0.5.*ones(size(u)) u], @(u,v)[u -0.5.*ones(size(u)) v], ...
          @(u,v)[ 0.5.*ones(size(u)) u v], @(u,v)[-0.5.*ones(size(u)) v u] };
for k = 1:6
    cube = cube + em.field.flux(FD3, faces{k}, [-0.5 0.5], [-0.5 0.5], 60, 60);
end
hemi = @(th,ph) 0.5.*[sin(th).*cos(ph) sin(th).*sin(ph) cos(th)];
capD = @(u,v) [u.*cos(v) -u.*sin(v) zeros(size(u))];
PhiHemi = em.field.flux(FD3, hemi, [0 pi/2], [0 2*pi], 60, 120, 'simpson') ...
        + em.field.flux(FD3, capD, [0 0.5], [0 2*pi], 60, 120, 'simpson');
fprintf('sphere              %.6e C\n', PhiSph);
fprintf('cube                %.6e C\n', cube);
fprintf('hemisphere and cap  %.6e C\n', PhiHemi);
em.test.assertClose(PhiSph, Q, 1e-4, 'Flux through the off-center sphere')
em.test.assertClose(cube, Q, 1e-4, 'Flux through the off-center cube')
em.test.assertClose(PhiHemi, Q, 1e-4, 'Flux through the off-center hemi-disk')

%% Independence of the charge position inside
for x0 = [0 0.15 0.3]
    sm = em.src.pointCharge(Q, [x0 0.1 0]);
    Pm = em.field.flux(@(r) em.field.D(sm, r), sphR, [0 pi], [0 2*pi], ...
        60, 120, 'simpson');
    fprintf('charge at x = %.2f   flux = %.6e C\n', x0, Pm);
    em.test.assertClose(Pm, Q, 1e-6, 'Flux at charge position')
end

%% Divergence against the density, an element count study
% A uniform ball of charge. The derivative step h must average over many
% elements, so h = 0.05 m is passed explicitly, and the element count is
% raised at fixed h.
rhov = 1e-6;
hdiv = 0.05;
for Nu = [10 20 40]
    ball = @(u,v,w)[w.*sin(u).*cos(v) w.*sin(u).*sin(v) w.*cos(u)];
    sb = em.src.volCharge(rhov, ball, [0 pi], [0 2*pi], [0 0.5], ...
        Nu, 2.*Nu, Nu, 'simpson');
    dv = em.op.div(@(r) em.field.D(sb, r), [0.1 0.05 0.1], hdiv);
    rel_error = abs(rhov - dv);
    if rel_error > 0.05
        error('ERROR: Relative error falls outside of tolerance.')
    end
    fprintf(['Nu = %2d   div D inside = %.4e   rho_v = %.4e   ' ...
        'rel err = %.2e\n'], Nu, dv, rhov, abs(dv - rhov)/rhov);
end
dfree = em.op.div(FD, [2 1 1], hdiv);
fprintf('div D in charge free space = %.3f  (consistent with zero)\n', dfree);
em.test.assertClose(dfree, 0, 1e-6, 'Divergence in charge free space')

%% Orientation is a convention you must own
sphFlip = @(ph,th) 0.5.*[sin(th).*cos(ph) sin(th).*sin(ph) cos(th)];
PhiFlip = em.field.flux(FD, sphFlip, [0 2*pi], [0 pi], 120, 60, 'simpson');
fprintf('flux with swapped parameter order = %.6e C  (sign flipped)\n', PhiFlip);
em.test.assertClose(PhiFlip, -Q, 1e-4, 'Flux of swapped parameter order')