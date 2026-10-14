function [V, info] = poisson2d(G, rho, epsr, mask, Vfix, opts)
%POISSON2D Solve div(eps grad(V)) = -rho on a uniform 2-D grid.
%   [V,INFO] = em.solve.poisson2d(G,rho,epsr,mask,Vfix,opts)
%   uses the five-point finite-difference stencil. epsr is relative
%   permittivity, so eps = eps0*epsr. rho may be a scalar, a Ny-by-Nx
%   array, or a function handle rho(X,Y). epsr may be a scalar or a
%   constant array matching the grid (spatially varying epsr is reserved
%   for the next lab).
%
%   MASK is true at fixed-potential nodes. VFIX is read only at those nodes.
%   Unmasked outer-edge nodes use homogeneous Neumann conditions, applied
%   by reflecting the missing neighbour (doubling the inward coefficient).
%
%   opts.method is 'jacobi', 'gs', 'sor', or 'direct'. Iterative methods
%   stop when the infinity-norm residual divided by max(norm(b,inf),1) is
%   no greater than opts.tol, or after opts.maxit sweeps. opts.omega is the
%   SOR relaxation factor.
%
%   INFO contains method, iter, time, resHist, converged, and finalResidual.
%   resHist starts with the initial residual (iteration zero).

if nargin < 6 || isempty(opts), opts = struct(); end
if ~isstruct(G) || ~all(isfield(G, {'X','Y','hx','hy'}))
    error('G must be created by em.solve.grid2d');
end
if ~isequal(size(G.X), size(G.Y)) || size(G.X,1) < 3 || size(G.X,2) < 3
    error('G.X and G.Y must be matching arrays at least 3-by-3');
end
sz = size(G.X);
if ~isequal(size(mask), sz) || ~islogical(mask)
    error('mask must be a logical array with the same size as G.X');
end
if ~isequal(size(Vfix), sz)
    error('Vfix must have the same size as G.X');
end
if any(~isfinite(Vfix(mask)))
    error('Vfix must be finite at all fixed-potential nodes');
end
if ~any(mask(:))
    error('At least one fixed-potential node is required to remove the constant nullspace');
end

if ~isfield(opts,'method') || isempty(opts.method), opts.method = 'sor'; end
if ~isfield(opts,'tol') || isempty(opts.tol), opts.tol = 1e-8; end
if ~isfield(opts,'maxit') || isempty(opts.maxit), opts.maxit = 10000; end
if ~isfield(opts,'omega') || isempty(opts.omega), opts.omega = 1.5; end
method = lower(char(opts.method));
if ~any(strcmp(method, {'jacobi','gs','sor','direct'}))
    error('opts.method must be jacobi, gs, sor, or direct');
end
if ~isscalar(opts.tol) || ~isfinite(opts.tol) || opts.tol <= 0
    error('opts.tol must be a positive finite scalar');
end
if ~isscalar(opts.maxit) || opts.maxit < 1 || opts.maxit ~= floor(opts.maxit)
    error('opts.maxit must be a positive integer');
end
if strcmp(method,'sor') && (~isscalar(opts.omega) || ~isfinite(opts.omega) || ...
        opts.omega <= 0 || opts.omega >= 2)
    error('SOR requires 0 < opts.omega < 2');
end
t0 = tic;

% Expand and validate the source density.
if isa(rho,'function_handle')
    rhoGrid = rho(G.X, G.Y);
elseif isscalar(rho)
    rhoGrid = repmat(rho, sz);
else
    rhoGrid = rho;
end
if isscalar(rhoGrid), rhoGrid = repmat(rhoGrid,sz); end
if ~isequal(size(rhoGrid),sz) || any(~isfinite(rhoGrid(:)))
    error('rho must be scalar, finite Ny-by-Nx data, or rho(X,Y) returning that size');
end

% Expand and validate relative permittivity.
if isscalar(epsr)
    epsGrid = repmat(epsr,sz);
else
    epsGrid = epsr;
end
if ~isequal(size(epsGrid),sz) || any(~isfinite(epsGrid(:))) || any(epsGrid(:) <= 0)
    error('epsr must be positive finite scalar or an array matching G.X');
end
% Lab 7 is the uniform-permittivity case. Reject spatially varying epsilon
% rather than silently applying a stencil that is not div(eps grad(V)).
if max(epsGrid(:)) - min(epsGrid(:)) > 100*eps(max(abs(epsGrid(:))))
    error('Lab 7 supports uniform epsr; piecewise permittivity is reserved for the next lab');
end

% Assemble A*Vfree=b. Each row represents -Laplacian(V)=rho/eps.
% Missing neighbours at unmasked outer edges are handled by homogeneous
% Neumann reflection, yielding twice the inward-neighbour coefficient.
freeLin = find(~mask);
nFree = numel(freeLin);
if nFree == 0
    V = Vfix;
    info = struct('method',method,'iter',0,'time',toc(t0),'resHist',0, ...
        'converged',true,'finalResidual',0);
    return
end
id = zeros(sz);
id(freeLin) = 1:nFree;

hx2i = 1/(G.hx^2);
hy2i = 1/(G.hy^2);
diagA = zeros(nFree,1);
b = zeros(nFree,1);
nbrId = zeros(nFree,4);
nbrC = zeros(nFree,4);

for k = 1:nFree
    lin = freeLin(k);
    [iy,ix] = ind2sub(sz,lin);
    b(k) = rhoGrid(lin)/(em.const.eps0()*epsGrid(lin));
    nn = 0;

    % x neighbours (columns).
    if ix == 1
        candidates = [iy, ix+1, 2*hx2i];
    elseif ix == sz(2)
        candidates = [iy, ix-1, 2*hx2i];
    else
        candidates = [iy, ix-1, hx2i; iy, ix+1, hx2i];
    end
    for q = 1:size(candidates,1)
        jy = candidates(q,1); jx = candidates(q,2); c = candidates(q,3);
        nn = nn+1; diagA(k) = diagA(k)+c;
        jlin = sub2ind(sz,jy,jx);
        if mask(jlin)
            b(k) = b(k) + c*Vfix(jlin);
        else
            nbrId(k,nn) = id(jlin); nbrC(k,nn) = c;
        end
    end

    % y neighbours (rows).
    if iy == 1
        candidates = [iy+1, ix, 2*hy2i];
    elseif iy == sz(1)
        candidates = [iy-1, ix, 2*hy2i];
    else
        candidates = [iy-1, ix, hy2i; iy+1, ix, hy2i];
    end
    for q = 1:size(candidates,1)
        jy = candidates(q,1); jx = candidates(q,2); c = candidates(q,3);
        nn = nn+1; diagA(k) = diagA(k)+c;
        jlin = sub2ind(sz,jy,jx);
        if mask(jlin)
            b(k) = b(k) + c*Vfix(jlin);
        else
            nbrId(k,nn) = id(jlin); nbrC(k,nn) = c;
        end
    end
end

% Sparse matrix for direct solution and residual evaluation.
I = zeros(nFree*5,1); J = zeros(nFree*5,1); S = zeros(nFree*5,1);
ptr = 0;
for k = 1:nFree
    ptr=ptr+1; I(ptr)=k; J(ptr)=k; S(ptr)=diagA(k);
    for q = 1:4
        j = nbrId(k,q);
        if j > 0
            ptr=ptr+1; I(ptr)=k; J(ptr)=j; S(ptr)=-nbrC(k,q);
        end
    end
end
A = sparse(I(1:ptr),J(1:ptr),S(1:ptr),nFree,nFree);
resScale = max(norm(b,inf),1);

V = Vfix;
x = zeros(nFree,1);
resHist = zeros(opts.maxit+1,1);
resHist(1) = norm(A*x-b,inf)/resScale;
converged = resHist(1) <= opts.tol;
iter = 0;

if strcmp(method,'direct')
    x = A\b;
    iter = 1;
    resHist = [resHist(1); norm(A*x-b,inf)/resScale];
    converged = resHist(end) <= max(opts.tol,1e-12);
else
    for it = 1:opts.maxit
        xOld = x;
        if strcmp(method,'jacobi')
            sums = zeros(nFree,1);
            for q = 1:4
                valid = nbrId(:,q) > 0;
                rows = find(valid);
                sums(rows) = sums(rows) + nbrC(rows,q).*xOld(nbrId(rows,q));
            end
            x = (b+sums)./diagA;
        else
            for k = 1:nFree
                ssum = 0;
                for q = 1:4
                    j = nbrId(k,q);
                    if j > 0
                        ssum = ssum + nbrC(k,q)*x(j);
                    end
                end
                xGS = (b(k)+ssum)/diagA(k);
                if strcmp(method,'sor')
                    x(k) = xOld(k) + opts.omega*(xGS-xOld(k));
                else
                    x(k) = xGS;
                end
            end
        end
        iter = it;
        resHist(it+1) = norm(A*x-b,inf)/resScale;
        if resHist(it+1) <= opts.tol
            converged = true;
            break
        end
    end
    resHist = resHist(1:iter+1);
end

V(freeLin) = x;
info = struct();
info.method = method;
info.iter = iter;
info.time = toc(t0);
info.resHist = resHist;
info.converged = converged;
info.finalResidual = resHist(end);
if ~converged && ~strcmp(method,'direct')
    warning('em:solve:poisson2d:NoConvergence', ...
        '%s did not reach tolerance %.3g after %d iterations (residual %.3g)', ...
        method, opts.tol, iter, info.finalResidual);
end
end
