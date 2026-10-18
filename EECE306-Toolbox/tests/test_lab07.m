function test_07 = test_lab07()
%TEST_LAB07 Lab 7 Test Function
%   Runs required tests for Lab 7

%% Setup
clear; close all;
e0 = em.const.eps0();

%% The manufactured solution
Vex = @(X,Y) sin(pi*X) .* sin(pi*Y);
rho = em.test.mms2d(Vex, 1);
G = em.solve.grid2d([0 1], [0 1], 81, 81);
mask = false(81,81); mask(1,:)=true; mask(end,:)=true; mask(:,1)=true; mask(:,end)=true;
Vfix = zeros(81,81);
[Vs, ~] = em.solve.poisson2d(G, rho, 1, mask, Vfix, struct('method','direct'));
err81 = max(abs(Vs(:) - Vex(G.X(:), G.Y(:))));
fprintf('MMS max error on the 81 x 81 grid = %.3e\n', err81);
if err81 > 1e-3
    error('ERROR: The max error of MMS exceeds 10^(-3) on a 101x101 grid.')
end

%% Convergence under grid refinement
Ns = [21 41 81 161]';
errs = zeros(size(Ns));
for j = 1:numel(Ns)
    n = Ns(j);
    Gn = em.solve.grid2d([0 1], [0 1], n, n);
    mk = false(n,n); mk(1,:)=true; mk(end,:)=true; mk(:,1)=true; mk(:,end)=true;
    Vn = em.solve.poisson2d(Gn, rho, 1, mk, zeros(n,n), struct('method','direct'));
    errs(j) = max(abs(Vn(:) - Vex(Gn.X(:), Gn.Y(:))));
    fprintf('N = %3d   max error = %.3e\n', n, errs(j));
end
p = em.test.convergence(@(N) errs(Ns == N), Ns);
fprintf('observed order p = %.2f  (expect near 2 for the five point stencil)\n', p);
em.test.assertClose(p, 2, 0.2, 'Observed order of convergence')

%% Three iterations and a direct solve
n = 41;
Gn = em.solve.grid2d([0 1], [0 1], n, n);
mk = false(n,n); mk(1,:)=true; mk(end,:)=true; mk(:,1)=true; mk(:,end)=true;
Vf = zeros(n,n);
methods = {'jacobi', 'gs', 'sor'};
maxits  = [3000 10000 10000];
figure; hold on
for k = 1:3
    o = struct('method', methods{k}, 'tol', 1e-6, 'maxit', maxits(k));
    if strcmp(methods{k}, 'sor'), o.omega = 2/(1+sin(pi/(n-1))); end
    [~, ik] = em.solve.poisson2d(Gn, rho, 1, mk, Vf, o);
    semilogy(ik.resHist)
    fprintf('%-7s iterations = %5d\n', methods{k}, ik.iter);
end
set(gca, 'YScale', 'log'); grid on; legend(methods)
xlabel('iteration'); ylabel('residual')
title('Three iterative methods on the same 41 x 41 problem')

%% The omega scan on two grids
figure; hold on
for n = [31 61]
    Gn = em.solve.grid2d([0 1], [0 1], n, n);
    mk = false(n,n); mk(1,:)=true; mk(end,:)=true; mk(:,1)=true; mk(:,end)=true;
    om = [1.0 1.4 1.7 1.85 1.90 1.95 1.99];
    it = zeros(size(om));
    for j = 1:numel(om)
        [~, ij] = em.solve.poisson2d(Gn, rho, 1, mk, zeros(n,n), ...
            struct('method','sor','tol',1e-4,'maxit',20000,'omega',om(j)));
        it(j) = ij.iter;
    end
    plot(om, it, 'o-')
    [~, jb] = min(it);
    fprintf('n = %3d   best omega about %.2f   theory 2/(1+sin(pi h)) = %.3f\n', ...
        n, om(jb), 2/(1+sin(pi/(n-1))));
end
grid on; legend('31 x 31', '61 x 61')
xlabel('omega'); ylabel('iterations to tolerance')
title('SOR is sharp, the best omega moves with the grid')

%% Timing, iterative against direct
Nt = [41 61 81];
td = zeros(size(Nt)); ts = zeros(size(Nt));
for k = 1:numel(Nt)
    n = Nt(k);
    Gn = em.solve.grid2d([0 1], [0 1], n, n);
    mk = false(n,n); mk(1,:)=true; mk(end,:)=true; mk(:,1)=true; mk(:,end)=true;
    tic; em.solve.poisson2d(Gn, rho, 1, mk, zeros(n,n), struct('method','direct')); td(k) = toc;
    tic; em.solve.poisson2d(Gn, rho, 1, mk, zeros(n,n), ...
        struct('method','sor','tol',1e-6,'maxit',40000,'omega',2/(1+sin(pi/(n-1))))); ts(k) = toc;
    fprintf('N = %3d   direct %.3f s   SOR %.3f s\n', n, td(k), ts(k));
end
figure; loglog(Nt, td, 'o-', Nt, ts, 's-'); grid on
xlabel('grid size N'); ylabel('wall time (s)')
legend('direct', 'SOR at optimal omega', 'location', 'northwest')
title('Cost against grid size')

%% A Laplace sanity pair
n = 41;
Gn = em.solve.grid2d([0 1], [0 1], n, n);
mk = false(n,n); mk(:,1) = true; mk(:,end) = true;
Vf = zeros(n,n); Vf(:,1) = 1;
VL = em.solve.poisson2d(Gn, 0, 1, mk, Vf, struct('method','direct'));
ramp = repmat(linspace(1, 0, n), n, 1);
fprintf('max deviation from the linear ramp = %.3e\n', max(abs(VL(:) - ramp(:))));
inner = VL(2:end-1, 2:end-1);
fprintf('interior range [%.4f, %.4f] inside boundary range [0, 1] = %d\n', ...
    min(inner(:)), max(inner(:)), min(inner(:)) >= 0 && max(inner(:)) <= 1);

%% Gauss's law closes the loop
n = 121;
Gn = em.solve.grid2d([-1 1], [-1 1], n, n);
mk = false(n,n); mk(1,:)=true; mk(end,:)=true; mk(:,1)=true; mk(:,end)=true;
rhoB = @(X,Y) 1e-8 * (abs(X) < 0.25 & abs(Y) < 0.25);
VB = em.solve.poisson2d(Gn, rhoB, 1, mk, zeros(n,n), struct('method','direct'));
Vh = @(r) interp2(Gn.X, Gn.Y, VB, r(:,1), r(:,2), 'linear');
Eh = @(r) -em.op.grad(Vh, r, Gn.hx/2);
edges = { @(t)[t -0.6+0*t 0*t], [0 -1 0]; @(t)[0.6+0*t t 0*t], [1 0 0]; ...
          @(t)[t 0.6+0*t 0*t], [0 1 0];  @(t)[-0.6+0*t t 0*t], [-1 0 0] };
flux = 0;
for k = 1:4
    [pos, w] = em.quad.line(edges{k,1}, [-0.6 0.6], 200);
    En = Eh(pos) * edges{k,2}';
    flux = flux + sum(En .* w);
end
Qenc = 1e-8 * 0.5 * 0.5;
fprintf('flux of E through the rectangle = %.4e\n', flux);
fprintf('enclosed charge over eps0       = %.4e\n', Qenc/e0);
fprintf('relative difference             = %.2e\n', abs(flux - Qenc/e0)/(Qenc/e0));
