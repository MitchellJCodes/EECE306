function [V, info] = poisson2d(G, rho, epsr, mask, Vfix, opts)
%POISSON2D Solve div(epsr grad(V)) = -rho/eps0 on a Cartesian grid.
% opts.method: jacobi|gs|sor|direct; opts.faceMean: harmonic (default)|arithmetic.
% Non-fixed outer edges have homogeneous Neumann conditions.

if nargin < 6 || isempty(opts), opts = struct(); end
if ~isfield(opts,'method'), opts.method = 'sor'; end
if ~isfield(opts,'tol'), opts.tol = 1e-8; end
if ~isfield(opts,'maxit'), opts.maxit = 20000; end
if ~isfield(opts,'omega'), opts.omega = 1.9; end
if ~isfield(opts,'faceMean'), opts.faceMean = 'harmonic'; end
if ~any(strcmpi(opts.faceMean, {'harmonic','arithmetic'}))
    error('poisson2d:badFaceMean','faceMean must be harmonic or arithmetic.');
end
if ~isscalar(opts.tol) || opts.tol <= 0 || ~isfinite(opts.tol)
    error('poisson2d:badTolerance','tol must be positive and finite.');
end
if ~isscalar(opts.maxit) || opts.maxit < 1 || opts.maxit ~= fix(opts.maxit)
    error('poisson2d:badMaxit','maxit must be a positive integer.');
end

[Ny,Nx] = size(G.X);
if Nx < 2 || Ny < 2, error('poisson2d:badGrid','Grid must be at least 2-by-2.'); end
R = expand(rho,Ny,Nx,G);
E = expand(epsr,Ny,Nx,G);
if ~isequal(size(R),[Ny Nx]) || ~isequal(size(E),[Ny Nx])
    error('poisson2d:badSize','rho and epsr must be scalar, matching arrays, or compatible handles.');
end
if any(~isfinite(E(:))) || any(E(:) <= 0)
    error('poisson2d:badPermittivity','epsr must be positive and finite.');
end
if isscalar(mask), mask = false(Ny,Nx); end
if isscalar(Vfix), Vfix = repmat(Vfix,Ny,Nx); end
if ~isequal(size(mask),[Ny Nx]) || ~isequal(size(Vfix),[Ny Nx])
    error('poisson2d:badSize','mask and Vfix must match the grid size.');
end
mask = logical(mask);
S = -R ./ em.const.eps0();

% Face coefficients are relative permittivity divided by spacing squared.
Ex = faceMean(E(:,1:end-1),E(:,2:end),opts.faceMean) / G.hx^2;
Ey = faceMean(E(1:end-1,:),E(2:end,:),opts.faceMean) / G.hy^2;

t0 = tic;
method = lower(opts.method);
switch method
    case 'direct'
        [V,A,b] = solveDirect(Ny,Nx,Ex,Ey,S,mask,Vfix);
        res = norm(A*V(:)-b,inf);
        info = struct('iter',1,'resHist',[],'time',toc(t0), ...
            'method','direct','finalResidual',res,'converged',isfinite(res));
    case {'jacobi','gs','sor'}
        [V,it,resHist,converged] = solveIterative(Ny,Nx,Ex,Ey,S,mask,Vfix,opts);
        info = struct('iter',it,'resHist',resHist,'time',toc(t0), ...
            'method',method,'finalResidual',resHist(end),'converged',converged);
    otherwise
        error('poisson2d:badMethod','Unknown method %s.',opts.method);
end
end

function [V,it,resHist,converged] = solveIterative(Ny,Nx,Ex,Ey,S,mask,Vfix,opts)
V = zeros(Ny,Nx); V(mask)=Vfix(mask);
jacobi = strcmpi(opts.method,'jacobi'); sor=strcmpi(opts.method,'sor');
resHist=zeros(opts.maxit,1); converged=false;
for it=1:opts.maxit
    Vold=V;
    for j=1:Nx
        for i=1:Ny
            if mask(i,j), continue; end
            num=S(i,j); den=0;
            if j>1, num=num-Ex(i,j-1)*getval(i,j-1,V,Vold,jacobi); den=den+Ex(i,j-1); end
            if j<Nx, num=num-Ex(i,j)*getval(i,j+1,V,Vold,jacobi); den=den+Ex(i,j); end
            if i>1, num=num-Ey(i-1,j)*getval(i-1,j,V,Vold,jacobi); den=den+Ey(i-1,j); end
            if i<Ny, num=num-Ey(i,j)*getval(i+1,j,V,Vold,jacobi); den=den+Ey(i,j); end
            vnew=num/(-den);
            if sor, V(i,j)=V(i,j)+opts.omega*(vnew-V(i,j)); else, V(i,j)=vnew; end
        end
    end
    r=norm(V(:)-Vold(:),inf); resHist(it)=r;
    if r<opts.tol, converged=true; break; end
end
resHist=resHist(1:it);
end

function v=getval(i,j,V,Vold,jacobi)
if jacobi, v=Vold(i,j); else, v=V(i,j); end
end

function [V,A,b] = solveDirect(Ny,Nx,Ex,Ey,S,mask,Vfix)
n=Ny*Nx; idx=@(i,j)(j-1)*Ny+i;
I=zeros(5*n,1); J=I; X=I; b=zeros(n,1); k=0;
for j=1:Nx
    for i=1:Ny
        q=idx(i,j);
        if mask(i,j)
            k=k+1; I(k)=q; J(k)=q; X(k)=1; b(q)=Vfix(i,j); continue;
        end
        den=0;
        if j>1
            a=Ex(i,j-1); den=den+a; k=k+1; I(k)=q; J(k)=idx(i,j-1); X(k)=a;
        end
        if j<Nx
            a=Ex(i,j); den=den+a; k=k+1; I(k)=q; J(k)=idx(i,j+1); X(k)=a;
        end
        if i>1
            a=Ey(i-1,j); den=den+a; k=k+1; I(k)=q; J(k)=idx(i-1,j); X(k)=a;
        end
        if i<Ny
            a=Ey(i,j); den=den+a; k=k+1; I(k)=q; J(k)=idx(i+1,j); X(k)=a;
        end
        k=k+1; I(k)=q; J(k)=q; X(k)=-den; b(q)=S(i,j);
    end
end
A=sparse(I(1:k),J(1:k),X(1:k),n,n);
V=reshape(A\b,Ny,Nx);
end

function F=faceMean(a,b,kind)
if strcmpi(kind,'arithmetic'), F=(a+b)/2;
else, F=2.*a.*b./(a+b); end
end

function M=expand(v,Ny,Nx,G)
if isa(v,'function_handle'), M=v(G.X,G.Y);
elseif isscalar(v), M=repmat(v,Ny,Nx);
else, M=v;
end
end
