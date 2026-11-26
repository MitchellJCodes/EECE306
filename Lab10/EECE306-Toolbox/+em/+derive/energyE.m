function W = energyE(G,V,epsr,opts)
%ENERGYE Stored electric energy per unit length in a 2-D cross-section.
% W = em.derive.energyE(G,V,epsr) uses harmonic face permittivity.
% Optional opts.faceMean selects 'harmonic' (default) or 'arithmetic'.
if nargin<4 || isempty(opts), opts=struct(); end
if ~isfield(opts,'faceMean'), opts.faceMean='harmonic'; end
[Ny,Nx]=size(G.X);
if ~isequal(size(V),[Ny Nx]), error('energyE:size','V must match G.'); end
if isscalar(epsr), E=repmat(epsr,Ny,Nx); else, E=epsr; end
if ~isequal(size(E),[Ny Nx]) || any(E(:)<=0) || any(~isfinite(E(:)))
    error('energyE:epsr','epsr must be positive and match G.');
end
Ex=faceMean(E(:,1:end-1),E(:,2:end),opts.faceMean);
Ey=faceMean(E(1:end-1,:),E(2:end,:),opts.faceMean);
wy=G.hy*ones(Ny,1); wy([1 end])=G.hy/2;
wx=G.hx*ones(1,Nx); wx([1 end])=G.hx/2;
Wx=sum(sum(Ex.*(diff(V,1,2).^2/G.hx^2).*repmat(wy,1,Nx-1)))*G.hx;
Wy=sum(sum(Ey.*(diff(V,1,1).^2/G.hy^2).*repmat(wx,Ny-1,1)))*G.hy;
W=0.5*em.const.eps0()*(Wx+Wy);
end
function F=faceMean(a,b,kind)
if strcmpi(kind,'arithmetic'), F=(a+b)/2;
elseif strcmpi(kind,'harmonic'), F=2.*a.*b./(a+b);
else, error('energyE:faceMean','faceMean must be harmonic or arithmetic.'); end
end
