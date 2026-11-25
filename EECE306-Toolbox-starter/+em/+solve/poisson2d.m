function [V, info] = poisson2d(G, rho, epsr, mask, Vfix, opts)
%POISSON2D  Solve div(eps grad V) = -rho on a 2D grid.  [Lab 7 skeleton]
%
%   Save this file as  +em/+solve/poisson2d.m  inside your EECE306-Toolbox.
%
%   The iteration loop, the convergence check, the boundary handling, the
%   direct sparse solve and the info bookkeeping are written for you. Fill
%   in the two places marked TODO, the stencil coefficients and the update
%   rule. Everything else, including the direct solver, runs off those two.
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
% Permitivity averaged at interfaces
if ~isfield(opts,'avgMethod'), opts.avgMethod = 'harmonic'; end

[Ny, Nx] = size(G.X);
hx2 = G.hx^2;
hy2 = G.hy^2;
e0  = em.const.eps0();

R = expandGrid(rho, Ny, Nx, G);
E = expandGrid(epsr, Ny, Nx, G);
% Checks permittivity is physically valid
if any(~isfinite(E(:))) || any(E(:) <= 0) 
    error('poisson2d:badEps', 'Relative permittivity must be positive and finite.');
end

if isscalar(mask), mask = false(Ny, Nx); end
if isscalar(Vfix), Vfix = repmat(Vfix, Ny, Nx); end

if ~any(mask(:))
    error('poisoon2d:noReference', 'At least one fixed-potential node is required.');
end

% Right-hand side of the stencil equation at every node, Laplacian(V) = S.
S = -R ./ e0;
% Placing the flag over how the permittivity is averaged between
% neighboring nodes
switch lower(opts.avgMethod)
    case 'harmonic'
        avg = @(e1,e2) 2.*e1.*e2./(e1+e2);

    case 'arithmetic'
        avg = @(e1,e2) (e1+e2)./2;
    otherwise
        error('poisson2d:badWverage','Unknown averaging method: %s', opts.avgMethod);
end
% ---- TODO 1 -----------------------------------------------------------
% Write the three stencil coefficients. Start from the second central
% difference in x and in y, add them, and set the sum equal to S(i,j).
% The equation at an interior node (i,j) must read
%
%   cx*(V(i,j-1) + V(i,j+1)) + cy*(V(i-1,j) + V(i+1,j)) + cc*V(i,j) = S(i,j)
%
% cx multiplies the two x neighbours, cy the two y neighbours, and cc is
% the coefficient of the centre node. This form is valid for the UNIFORM
% permittivity used in this lab. Lab 8 upgrades the stencil so each of the
% four faces carries its own averaged eps, which a dielectric interface
% needs. Keep that in mind when you choose how to structure the code.
% cx = 1 / hx2;   % <-- coefficient x
% cy = 1 / hy2;   % <-- coefficient y
% cc = -2 * (cx + cy);   % <-- coefficient c
% -----------------------------------------------------------------------

t0 = tic;
switch lower(opts.method)
    case 'direct'
        V = solveDirect(Ny, Nx, E, hx2, hy2, avg, S, mask, Vfix);
        info = struct('iter', 1, 'resHist', [], 'time', toc(t0), 'method', 'direct');
    case {'jacobi', 'gs', 'sor'}
        [V, it, resHist] = solveIterative(Ny, Nx, E, hx2, hy2, avg, S, mask, Vfix, opts);
        info = struct('iter', it, 'resHist', resHist, 'time', toc(t0), 'method', lower(opts.method));
    otherwise
        error('poisson2d:badMethod', 'Unknown method %s.', opts.method);
end
end

% ---- Jacobi, Gauss-Seidel, SOR --------------------------------------------
function [V, it, resHist] = solveIterative(Ny, Nx, E, hx2, hy2, avg, S, mask, Vfix, opts)
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
            
            % calculate permittivityp at each face using
            % harmonic/arithmetic averaging
            [ce, cw, cn, cs, cc] = faceCoeffs(E, i, j, iu, id, jl, jr, hx2, hy2, avg);
            % Jacobi reads every neighbour from the previous sweep.
            % Gauss-Seidel and SOR read the newest value available.
            if jacobi
                vu = Vold(iu,j); vd = Vold(id,j); vl = Vold(i,jl); vr = Vold(i,jr);
            else
                vu = V(iu,j);    vd = V(id,j);    vl = V(i,jl);    vr = V(i,jr);
            end

            % ---- TODO 2 -----------------------------------------------
            % Solve the stencil equation of TODO 1 for the centre value.
            % Use vl, vr (x neighbours), vu, vd (y neighbours), S(i,j),
            % and the coefficients cx, cy, cc.
            % vnew = (S(i,j) - cx*(vl + vr) - cy*(vu + vd)) / cc;
            vnew = (S(i,j) - ce*vr - cw*vl - cn*vu - cs*vd) / cc;
            % -----------------------------------------------------------

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
% function V = solveDirect(Ny, Nx, cx, cy, cc, S, mask, Vfix)
% The same stencil, assembled once as a sparse matrix and solved with
% backslash. Duplicate entries in sparse() add up, which is exactly what
% the reflected neighbour at a Neumann edge needs.
% n   = Ny * Nx;
% idx = @(i, j) (j - 1) * Ny + i;
% I = zeros(5*n, 1); J = I; X = I; b = zeros(n, 1); k = 0;
% for j = 1:Nx
%     for i = 1:Ny
%         p = idx(i, j);
%         if mask(i,j)
%             k = k + 1; I(k) = p; J(k) = p; X(k) = 1; b(p) = Vfix(i,j);
%             continue;
%         end
%         [iu, id, jl, jr] = neigh(i, j, Ny, Nx);
%         k = k + 1; I(k) = p; J(k) = p;          X(k) = cc;
%         k = k + 1; I(k) = p; J(k) = idx(i, jr); X(k) = cx;
%         k = k + 1; I(k) = p; J(k) = idx(i, jl); X(k) = cx;
%         k = k + 1; I(k) = p; J(k) = idx(iu, j); X(k) = cy;
%         k = k + 1; I(k) = p; J(k) = idx(id, j); X(k) = cy;
%         b(p) = S(i,j);
%     end
% end
% A = sparse(I(1:k), J(1:k), X(1:k), n, n);
% V = reshape(A \ b, Ny, Nx);
% end
function V = solveDirect(Ny, Nx, E, hx2, hy2, avg, S, mask, Vfix)
% Assemble a sparse matrix for the
% dielectric-aware five-point stencil.

n = Ny * Nx;

idx = @(i,j) (j-1)*Ny + i;

I = zeros(5*n,1);
J = zeros(5*n,1);
X = zeros(5*n,1);

b = zeros(n,1);

k = 0;

for j = 1:Nx
    for i = 1:Ny

        p = idx(i,j);

        % Apply fixed-potential boundary conditions
        if mask(i,j)

            k = k + 1;

            I(k) = p;
            J(k) = p;
            X(k) = 1;

            b(p) = Vfix(i,j);

            continue;
        end
        % Find neighboring grid indices
        [iu, id, jl, jr] = neigh(i, j, Ny, Nx);

        % Calculate dielectric-aware stencil coefficients
        [ce, cw, cn, cs, cc] = faceCoeffs( ...
            E, i, j, iu, id, jl, jr, hx2, hy2, avg);

        % Center coefficient
        k = k + 1;
        I(k) = p;
        J(k) = p;
        X(k) = cc;

        % East neighbor
        k = k + 1;
        I(k) = p;
        J(k) = idx(i,jr);
        X(k) = ce;

        % West neighbor
        k = k + 1;
        I(k) = p;
        J(k) = idx(i,jl);
        X(k) = cw;

        % North neighbor
        k = k + 1;
        I(k) = p;
        J(k) = idx(iu,j);
        X(k) = cn;

        % South neighbor
        k = k + 1;
        I(k) = p;
        J(k) = idx(id,j);
        X(k) = cs;

        % Source term
        b(p) = S(i,j);

    end
end

% Assemble sparse matrix
A = sparse(I(1:k), J(1:k), X(1:k), n, n);

% Solve the linear system
V = reshape(A\b, Ny, Nx);

end
% ----------- Face-dependent stencil coefficient --------------------------
function [ce, cw, cn, cs, cc] = faceCoeffs( ...
    E, i, j, iu, id, jl, jr, hx2, hy2, avg)

% Calculate averaged relative permittivity
% between the current node and its neighbors.

epsE = avg(E(i,j), E(i,jr));
epsW = avg(E(i,j), E(i,jl));

epsN = avg(E(i,j), E(iu,j));
epsS = avg(E(i,j), E(id,j));

% Calculate four face coefficients
ce = epsE / hx2;
cw = epsW / hx2;

cn = epsN / hy2;
cs = epsS / hy2;

% Center coefficient
cc = -(ce + cw + cn + cs);

end

% ---- helpers ---------------------------------------------------------------
function [iu, id, jl, jr] = neigh(i, j, Ny, Nx)

iu = min(i+1, Ny);
id = max(i-1, 1);
jr = min(j+1, Nx);
jl = max(j-1, 1);

if i == 1
    id = 2;
elseif i == Ny
    iu = Ny-1;
end

if j == 1
    jl = 2;
elseif j == Nx
    jr = Nx-1;
end

end

function M = expandGrid(v, Ny, Nx, G)

if isa(v, 'function_handle')
    M = v(G.X, G.Y);
elseif isscalar(v)
    M = repmat(v, Ny, Nx);
else
    M = v;
end

end