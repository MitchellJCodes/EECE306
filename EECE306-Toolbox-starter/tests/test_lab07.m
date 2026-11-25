function test_07 = test_lab07()
%TEST_LAB07 Verify the Lab 7 Poisson solver (seven required checks).
% Run from the toolbox root: addpath(pwd); test_lab07()
% Assumes grid arrays G.X/G.Y and potential V use the same indexing.
% Requires: grid2d, poisson2d, and the Lab 4 convergence helper.
%
% NOTE: Test 7 uses a direct numerical contour flux calculation so the
% result is independent of the exact calling convention of em.field.flux.
% To satisfy the handout fully, replace the indicated calculation with
% em.field.flux using the signature from your Lab 5 implementation.

ok = false;
eps0 = em.const.eps0();
Vexact = @(x,y) sin(pi*x).*sin(pi*y);
opts = struct('method','direct','tol',1e-9,'maxit',20000,'omega',1.8);

%% 1. Manufactured solution on a 101 x 101 grid
G = em.solve.grid2d([0 1],[0 1],101,101);
[mask,Vfix] = zeroBoundary(G);
rho = 2*pi^2*eps0*Vexact(G.X,G.Y);
[V,~] = em.solve.poisson2d(G,rho,1,mask,Vfix,opts);
err = max(abs(V(:)-Vexact(G.X(:),G.Y(:))));
assert(err < 1e-3, 'Test 1: MMS maximum error is too large');
fprintf('Test 1 passed: MMS error = %.3e\n',err);

%% 2. Second-order convergence of the five-point stencil
Ns = [11 21 41 81];
errors = zeros(size(Ns));

for k = 1:numel(Ns)
    g = em.solve.grid2d([0 1],[0 1],Ns(k),Ns(k));
    [m,vf] = zeroBoundary(g);
    rr = 2*pi^2*eps0*Vexact(g.X,g.Y);
    u = em.solve.poisson2d(g,rr,1,m,vf,opts);
    errors(k) = max(abs(u(:)-Vexact(g.X(:),g.Y(:))));
end

h = 1./(Ns-1);
p = polyfit(log(h),log(errors),1);
order = p(1);

assert(order > 1.8 && order < 2.2, ...
    'Test 2: order not near two');

fprintf('Test 2 passed: convergence order = %.3f\n',order);

%% 3 and 4. Four methods agree; iteration ordering
% Moderate grid so Jacobi finishes in a reasonable time.

g = em.solve.grid2d([0 1],[0 1],21,21);
[m,vf] = zeroBoundary(g);
rr = 2*pi^2*eps0*Vexact(g.X,g.Y);

methods = {'jacobi','gs','sor','direct'};
U = cell(1,4);
info = cell(1,4);

for k = 1:4
    o = struct('method',methods{k},'tol',1e-10, ...
               'maxit',30000,'omega',1.7);

    [U{k},info{k}] = em.solve.poisson2d(g,rr,1,m,vf,o);
end

% Compare iterative methods against the direct solution.
for k = 1:3
    delta = max(abs(U{k}(:)-U{4}(:)));

    assert(delta < 1e-6, ...
        'Test 3: solver methods disagree');
end

fprintf('Test 3 passed: Jacobi, GS, SOR, direct agree\n');

assert(info{3}.iter < info{2}.iter && ...
       info{2}.iter < info{1}.iter, ...
       'Test 4: expected SOR < GS < Jacobi iteration counts');

fprintf('Test 4 passed: iterations SOR=%d, GS=%d, Jacobi=%d\n', ...
        info{3}.iter,info{2}.iter,info{1}.iter);

%% 5. Linear Laplace solution with mixed boundary conditions
% V(x,y) = 1-x
% Dirichlet on left/right; Neumann on top/bottom.

g = em.solve.grid2d([0 1],[0 1],41,41);

m = false(size(g.X));
vf = zeros(size(g.X));

xmin = min(g.X(:));
xmax = max(g.X(:));

m(abs(g.X-xmin)<1e-12 | abs(g.X-xmax)<1e-12) = true;
vf(abs(g.X-xmin)<1e-12) = 1;

o = struct('method','direct','tol',1e-10, ...
           'maxit',30000,'omega',1.8);

Ulinear = em.solve.poisson2d(g,0,1,m,vf,o);

linearErr = max(abs(Ulinear(:)-(1-g.X(:))));

assert(linearErr < 1e-3, ...
    'Test 5: linear Laplace profile failed');

fprintf('Test 5 passed: mixed BC linear profile error %.3e\n', ...
        linearErr);

%% 6. Discrete maximum principle (source-free interior)

interior = Ulinear(2:end-1,2:end-1);

assert(all(interior(:) >= min(Ulinear(m))-1e-6) && ...
       all(interior(:) <= max(Ulinear(m))+1e-6), ...
       'Test 6: interior potential outside boundary extrema');

fprintf(['Test 6 passed: no source-free interior extrema ', ...
         'beyond BC bounds\n']);

%% 7. Gauss's law for a field derived from solved potential
% Reuse the manufactured solution.
% Evaluate D = -eps0 grad(V) using central differences.
% Compare flux through a rectangle with enclosed charge.
%
% This numerical flux calculation checks Gauss's law but
% does not directly call em.field.flux.

G = em.solve.grid2d([0 1],[0 1],101,101);
[m,vf] = zeroBoundary(G);

rr = 2*pi^2*eps0*Vexact(G.X,G.Y);

U = em.solve.poisson2d(G,rr,1,m,vf,opts);

[ux,uy] = numericGradient(U,G);

Dx = -eps0*ux;
Dy = -eps0*uy;

% Integrate D.n over a rectangle in [0.2,0.8].
ix = find(abs(G.X(1,:)-0.2)<1e-10,1);

if isempty(ix)
    [flux,Q] = fluxNdgrid(G,Dx,Dy,rr);
else
    [flux,Q] = fluxMeshgrid(G,Dx,Dy,rr);
end

rel = abs(flux-Q)/max(abs(Q),eps);

assert(rel < 0.06, ...
    'Test 7: Gauss law flux differs from enclosed charge');

fprintf('Test 7 passed: Gauss law relative error = %.3e\n',rel);

fprintf('All seven Lab 7 tests passed!\n');

ok = true;

end

%% Helper Function 1: Boundary Conditions

function [mask,Vfix] = zeroBoundary(G)

mask = false(size(G.X));

mask([1 end],:) = true;
mask(:,[1 end]) = true;

Vfix = zeros(size(G.X));

end

%% Helper Function 2: Numerical Gradient

function [ux,uy] = numericGradient(U,G)

if abs(G.X(1,2)-G.X(1,1)) > 0
    [ux,uy] = gradient(U,G.hx,G.hy);
else
    [uy,ux] = gradient(U,G.hy,G.hx);
end

end

%% Helper Function 3: Flux for Meshgrid

function [flux,Q] = fluxMeshgrid(G,Dx,Dy,rho)

i1 = 21;
i2 = 81;
j1 = 21;
j2 = 81;

flux = G.hy*sum(Dx(j1:j2,i2)-Dx(j1:j2,i1)) + ...
       G.hx*sum(Dy(j2,i1:i2)-Dy(j1,i1:i2));

Q = G.hx*G.hy*sum(sum(rho(j1:j2,i1:i2)));

end

%% Helper Function 4: Flux for Ndgrid

function [flux,Q] = fluxNdgrid(G,Dx,Dy,rho)

i1 = 21;
i2 = 81;
j1 = 21;
j2 = 81;

flux = G.hy*sum(Dx(i2,j1:j2)-Dx(i1,j1:j2)) + ...
       G.hx*sum(Dy(i1:i2,j2)-Dy(i1:i2,j1));

Q = G.hx*G.hy*sum(sum(rho(i1:i2,j1:j2)));

end
