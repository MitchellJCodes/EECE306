function [V, info] = poisson2d(G, rho, epsr, mask, Vfix, opts)
%POISSON2D  Solve div(eps grad V) = -rho on a 2D grid.  [Lab 7 skeleton]
%
%   G      grid struct from em.solve.grid2d
%   rho    volume charge density, array, scalar, or handle @(X,Y)
%   epsr   relative permittivity, array or scalar (scalar in this lab)
%   mask   logical, true where the potential is fixed (Dirichlet)
%   Vfix   values at the fixed nodes, read only where mask is true
%   opts   .method jacobi|gs|sor|direct, .tol, .maxit, .omega
%
%   Outer boundary nodes that are not fixed take the homogeneous Neumann
%   condition dV/dn = 0, implemented by reflection in neigh() below.
%   info.resHist is the residual at every iteration, info.iter and
%   info.time record the cost.

if nargin < 6, opts = struct; end
if ~isfield(opts,'method'), opts.method = 'sor';  end
if ~isfield(opts,'tol'),    opts.tol    = 1e-8;   end
if ~isfield(opts,'maxit'),  opts.maxit  = 20000;  end
if ~isfield(opts,'omega'),  opts.omega  = 1.9;    end

[Ny, Nx] = size(G.X);
hx2 = G.hx^2;
hy2 = G.hy^2;
e0  = em.const.eps0();

R = expand(rho,  Ny, Nx, G);
E = expand(epsr, Ny, Nx, G);
if isscalar(mask), mask = false(Ny, Nx); end
if isscalar(Vfix), Vfix = repmat(Vfix, Ny, Nx); end

% Right-hand side of the stencil equation at every node, Laplacian(V) = S.
S = -R ./ (e0 * E);

cx = 1 / hx2;
cy = 1 / hy2;
cc = -2 * cx - 2 * cy;

if cc == 0
    error('poisson2d:todo', 'Fill in TODO 1, the stencil coefficients, before running the solver.');
end

t0 = tic;
switch lower(opts.method)
    case 'direct'
        V = solveDirect(Ny, Nx, cx, cy, cc, S, mask, Vfix);
        info = struct('iter', 1, 'resHist', [], 'time', toc(t0), 'method', 'direct');
    case {'jacobi', 'gs', 'sor'}
        [V, it, resHist] = solveIterative(Ny, Nx, cx, cy, cc, S, mask, Vfix, opts);
        info = struct('iter', it, 'resHist', resHist, 'time', toc(t0), 'method', lower(opts.method));
    otherwise
        error('poisson2d:badMethod', 'Unknown method %s.', opts.method);
end
end

% ---- Jacobi, Gauss-Seidel, SOR --------------------------------------------
function [V, it, resHist] = solveIterative(Ny, Nx, cx, cy, cc, S, mask, Vfix, opts)
V = zeros(Ny, Nx);
V(mask) = Vfix(mask);
jacobi = strcmpi(opts.method, 'jacobi');
sor    = strcmpi(opts.method, 'sor');
w      = opts.omega;
resHist = zeros(opts.maxit, 1);

for it = 1:opts.maxit
    Vold = V;
    for j = 1:Nx
        for i = 1:Ny
            if mask(i,j), continue; end
            [iu, id, jl, jr] = neigh(i, j, Ny, Nx);

            % Jacobi reads every neighbour from the previous sweep.
            % Gauss-Seidel and SOR read the newest value available.
            if jacobi
                vu = Vold(iu,j); vd = Vold(id,j); vl = Vold(i,jl); vr = Vold(i,jr);
            else
                vu = V(iu,j);    vd = V(id,j);    vl = V(i,jl);    vr = V(i,jr);
            end

            vnew = (S(i,j) - cx*(vl + vr) - cy*(vu + vd)) / cc;

            if sor
                V(i,j) = V(i,j) + w * (vnew - V(i,j));
            else
                V(i,j) = vnew;
            end
        end
    end
    r = norm(V(:) - Vold(:), inf);
    resHist(it) = r;
    if r < opts.tol, break; end
end
resHist = resHist(1:it);
end

% ---- direct sparse solve ---------------------------------------------------
function V = solveDirect(Ny, Nx, cx, cy, cc, S, mask, Vfix)
% The same stencil, assembled once as a sparse matrix and solved with
% backslash. Duplicate entries in sparse() add up, which is exactly what
% the reflected neighbour at a Neumann edge needs.
n   = Ny * Nx;
idx = @(i, j) (j - 1) * Ny + i;
I = zeros(5*n, 1); J = I; X = I; b = zeros(n, 1); k = 0;
for j = 1:Nx
    for i = 1:Ny
        p = idx(i, j);
        if mask(i,j)
            k = k + 1; I(k) = p; J(k) = p; X(k) = 1; b(p) = Vfix(i,j);
            continue;
        end
        [iu, id, jl, jr] = neigh(i, j, Ny, Nx);
        k = k + 1; I(k) = p; J(k) = p;          X(k) = cc;
        k = k + 1; I(k) = p; J(k) = idx(i, jr); X(k) = cx;
        k = k + 1; I(k) = p; J(k) = idx(i, jl); X(k) = cx;
        k = k + 1; I(k) = p; J(k) = idx(iu, j); X(k) = cy;
        k = k + 1; I(k) = p; J(k) = idx(id, j); X(k) = cy;
        b(p) = S(i,j);
    end
end
A = sparse(I(1:k), J(1:k), X(1:k), n, n);
V = reshape(A \ b, Ny, Nx);
end

% ---- helpers ---------------------------------------------------------------
function [iu, id, jl, jr] = neigh(i, j, Ny, Nx)
% Homogeneous Neumann at an outer edge, implemented by reflection.
iu = min(i+1, Ny);  id = max(i-1, 1);
jr = min(j+1, Nx);  jl = max(j-1, 1);
end

function M = expand(v, Ny, Nx, G)
if isa(v, 'function_handle'), M = v(G.X, G.Y);
elseif isscalar(v),           M = repmat(v, Ny, Nx);
else,                         M = v;
end
end
