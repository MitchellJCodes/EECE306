function test_07 = test_lab07()
%TEST_LAB07 Verification tests for the 2-D Poisson solver.

Vexact = @(X,Y) sin(pi*X).*sin(pi*Y);
rhoMMS = em.test.mms2d(Vexact,1);

%% Pointwise Laplacian operator
fLap = @(r) sin(pi*r(:,1)).*sin(pi*r(:,2));
rLap = [0.23 0.41 0; 0.55 0.72 0; 0.81 0.36 0];
lapNum = em.op.laplacian(fLap,rLap,1e-3);
lapRef = -2*pi^2*fLap(rLap);
assert(max(abs(lapNum-lapRef)) < 1e-4, ...
    'Pointwise Laplacian does not match the analytic Laplacian');

%% 1. Manufactured solution on a 101-by-101 grid
G = em.solve.grid2d([0 1],[0 1],101,101);
mask = false(size(G.X));
mask([1 end],:) = true; mask(:,[1 end]) = true;
Vfix = zeros(size(G.X));
[V,info] = em.solve.poisson2d(G,rhoMMS,1,mask,Vfix, ...
    struct('method','direct','tol',1e-10));
err = max(abs(V(:)-Vexact(G.X(:),G.Y(:))));
assert(err < 1e-3, sprintf('MMS max error %.3g exceeds 1e-3',err));
assert(info.finalResidual < 1e-8,'Direct solve residual is unexpectedly large');

%% 2. Observed second-order convergence
Ns = [21 41 81];
errs = zeros(size(Ns));
for k = 1:numel(Ns)
    n = Ns(k);
    Gk = em.solve.grid2d([0 1],[0 1],n,n);
    mk = false(size(Gk.X)); mk([1 end],:) = true; mk(:,[1 end]) = true;
    Vk = em.solve.poisson2d(Gk,rhoMMS,1,mk,zeros(size(Gk.X)), ...
        struct('method','direct'));
    errs(k) = max(abs(Vk(:)-Vexact(Gk.X(:),Gk.Y(:))));
end
[p,hfig] = em.test.convergence(@(n) errs(Ns == n),Ns, ...
    struct('title','Lab 7 manufactured-solution convergence'));
assert(p > 1.8 && p < 2.2, ...
    sprintf('Observed convergence order %.3f is outside [1.8,2.2]',p));
close(hfig);

%% 3. Jacobi, Gauss-Seidel, SOR, and direct agree
n = 21;
G = em.solve.grid2d([0 1],[0 1],n,n);
mask = false(size(G.X)); mask([1 end],:) = true; mask(:,[1 end]) = true;
Vfix = zeros(size(G.X));
common = struct('tol',1e-10,'maxit',12000);
optsDirect = struct('method','direct','tol',1e-10);
optsJacobi = common; optsJacobi.method = 'jacobi'; optsJacobi.omega = 1.5;
optsGS = common; optsGS.method = 'gs'; optsGS.omega = 1.5;
omegaCompare = 2 / (1 + sin(pi / (n - 1)));
optsSOR = common; optsSOR.method = 'sor'; optsSOR.omega = omegaCompare;
[Vd, id] = em.solve.poisson2d(G,rhoMMS,1,mask,Vfix,optsDirect);
[Vj, ij] = em.solve.poisson2d(G,rhoMMS,1,mask,Vfix,optsJacobi);
[Vg, ig] = em.solve.poisson2d(G,rhoMMS,1,mask,Vfix,optsGS);
[Vs, isor] = em.solve.poisson2d(G,rhoMMS,1,mask,Vfix,optsSOR);
assert(max(abs(Vj(:)-Vd(:))) < 2e-7,'Jacobi disagrees with direct solution');
assert(max(abs(Vg(:)-Vd(:))) < 2e-7,'Gauss-Seidel disagrees with direct solution');
assert(max(abs(Vs(:)-Vd(:))) < 2e-7,'SOR disagrees with direct solution');
assert(id.converged && ij.converged && ig.converged && isor.converged, ...
    'One or more solver methods did not converge');

%% 4. SOR is faster in iterations than GS, which is faster than Jacobi
n = 31;
G = em.solve.grid2d([0 1],[0 1],n,n);
mask = false(size(G.X)); mask([1 end],:) = true; mask(:,[1 end]) = true;
Vfix = zeros(size(G.X));
opt = 2/(1+sin(pi/(n-1)));
[~,ij] = em.solve.poisson2d(G,rhoMMS,1,mask,Vfix, ...
    struct('method','jacobi','tol',1e-6,'maxit',12000));
[~,ig] = em.solve.poisson2d(G,rhoMMS,1,mask,Vfix, ...
    struct('method','gs','tol',1e-6,'maxit',12000));
[~,isor] = em.solve.poisson2d(G,rhoMMS,1,mask,Vfix, ...
    struct('method','sor','tol',1e-6,'maxit',12000,'omega',opt));
assert(isor.iter < ig.iter && ig.iter < ij.iter, ...
    sprintf('Expected SOR < GS < Jacobi iterations; got %d, %d, %d', ...
    isor.iter,ig.iter,ij.iter));

%% 5. Laplace with Dirichlet on left/right and Neumann on top/bottom
G = em.solve.grid2d([0 1],[0 1],31,25);
mask = false(size(G.X)); mask(:,1)=true; mask(:,end)=true;
Vfix = zeros(size(G.X)); Vfix(:,1)=1; Vfix(:,end)=0;
VL = em.solve.poisson2d(G,0,1,mask,Vfix,struct('method','direct'));
ramp = repmat(linspace(1,0,G.Nx),G.Ny,1);
assert(max(abs(VL(:)-ramp(:))) < 1e-10, ...
    'Mixed-boundary Laplace problem did not produce a linear profile');

%% 6. A source-free solution has no interior extrema
U = VL(2:end-1,2:end-1);
N = cat(3,VL(1:end-2,2:end-1),VL(3:end,2:end-1), ...
    VL(2:end-1,1:end-2),VL(2:end-1,3:end));
localMin = min(N,[],3); localMax = max(N,[],3);
assert(all(U(:) >= localMin(:)-1e-12) && all(U(:) <= localMax(:)+1e-12), ...
    'Source-free potential contains an interior local extremum');

%% 7. Gauss law for a potential obtained from the solver
n = 101;
G = em.solve.grid2d([-1 1],[-1 1],n,n);
mask = false(size(G.X)); mask([1 end],:)=true; mask(:,[1 end])=true;
rhoBlock = @(X,Y) 1e-8 .* (abs(X)<0.25 & abs(Y)<0.25);
rhoGrid = rhoBlock(G.X,G.Y);
VB = em.solve.poisson2d(G,rhoBlock,1,mask,zeros(size(G.X)), ...
    struct('method','direct'));
Vh = @(r) interp2(G.X,G.Y,VB,r(:,1),r(:,2),'linear');
Eh = @(r) -em.op.grad(Vh,r,G.hx/2);
Dh = @(r) em.const.eps0()*Eh(r);
% Four oriented vertical faces form the side of a unit-height rectangular
% Gaussian box. The top/bottom flux is zero because this 2-D solution has
% no z dependence and Ez=0.
faces = { ...
    @(u,v)[u, -0.6+0*u, v], ... % outward -y
    @(u,v)[0.6+0*u, u, v], ... % outward +x
    @(u,v)[-u, 0.6+0*u, v], ... % outward +y
    @(u,v)[-0.6+0*u, -u, v]};  % outward -x
fluxTotal = 0;
for k = 1:numel(faces)
    fluxTotal = fluxTotal + em.field.flux(Dh,faces{k},[-0.6 0.6],[0 1],120,2);
end
Qenc = sum(rhoGrid(:))*G.hx*G.hy; % unit height, grid-aligned source block
fluxExpected = Qenc;
relFluxErr = abs(fluxTotal-fluxExpected)/max(abs(fluxExpected),eps);
assert(relFluxErr < 0.15, ...
    sprintf('Gauss-law relative error %.3g exceeds 15%%',relFluxErr));

fprintf('  Lab 7 MMS error %.3e; observed order %.3f; Gauss-law relative error %.3f\n', ...
    err,p,relFluxErr);
disp('  test_lab07 checks complete');
end

