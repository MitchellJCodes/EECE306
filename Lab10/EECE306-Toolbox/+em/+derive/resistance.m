function R = resistance(G,sigma,mask,Vfix,opts)
%RESISTANCE Compute resistance per unit length from total conduction current.
% Integrates J=-sigma grad(V) across a vertical cut through the grid.
if nargin<5 || isempty(opts), opts=struct(); end
[Ny,Nx]=size(G.X);
if isscalar(sigma), S=repmat(sigma,Ny,Nx); else, S=sigma; end
if ~isequal(size(S),[Ny Nx]) || any(S(:)<=0) || any(~isfinite(S(:)))
    error('resistance:sigma','sigma must be positive, finite, and match G.');
end
if isscalar(mask), mask=false(Ny,Nx); end
if isscalar(Vfix), Vfix=repmat(Vfix,Ny,Nx); end
vals=Vfix(logical(mask));
if numel(vals)<2 || max(vals)<=min(vals)
    error('resistance:conductors','mask must specify two distinct fixed potentials.');
end
[V,~]=em.solve.poisson2d(G,0,S,mask,Vfix,opts);
j=max(1,min(Nx-1,round(Nx/2)));
sigmaFace=2.*S(:,j).*S(:,j+1)./(S(:,j)+S(:,j+1));
Jx=-sigmaFace.*(V(:,j+1)-V(:,j))/G.hx;
weights=G.hy*ones(Ny,1); weights([1 end])=G.hy/2;
I=abs(sum(Jx.*weights));
if I<=0 || ~isfinite(I), error('resistance:current','Computed current is zero or invalid.'); end
V0=max(vals)-min(vals);
R=V0/I;
end
