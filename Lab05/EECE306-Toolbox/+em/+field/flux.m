function Phi = flux(F, surface, uinterval, vinterval, Nu, Nv)
%FLUX Compute the flux of a vector field through a parametrized surface.
%
%   Phi = em.field.flux(F,surface,uinterval,vinterval,Nu,Nv)
%
%   F          - vector field function handle
%   surface    - parametrized surface r(u,v)
%   uinterval  - [u0 u1]
%   vinterval  - [v0 v1]
%   Nu, Nv     - number of quadrature points
%
%   Phi        - signed flux through the surface

% Parameter-space quadrature
[u,wu] = em.quad.nodes(uinterval(1),uinterval(2),Nu,'midpoint');
[v,wv] = em.quad.nodes(vinterval(1),vinterval(2),Nv,'midpoint');

[U,V] = ndgrid(u,v);

U = U(:);
V = V(:);

[WU,WV] = ndgrid(wu,wv);
wt = WU(:).*WV(:);

% Surface positions
r = surface(U,V);
r = reshape(r,[],3);

% Numerical derivative steps
hu = eps^(1/5) .* (1 + abs(U));
hv = eps^(1/5) .* (1 + abs(V));

% Five-point derivative with respect to u
r_um2 = surface(U - 2*hu,V);
r_um1 = surface(U - hu,V);
r_up1 = surface(U + hu,V);
r_up2 = surface(U + 2*hu,V);

% Five-point derivative with respect to v
r_vm2 = surface(U,V - 2*hv);
r_vm1 = surface(U,V - hv);
r_vp1 = surface(U,V + hv);
r_vp2 = surface(U,V + 2*hv);

% Reshape derivative results
r_um2 = reshape(r_um2,[],3);
r_um1 = reshape(r_um1,[],3);
r_up1 = reshape(r_up1,[],3);
r_up2 = reshape(r_up2,[],3);

r_vm2 = reshape(r_vm2,[],3);
r_vm1 = reshape(r_vm1,[],3);
r_vp1 = reshape(r_vp1,[],3);
r_vp2 = reshape(r_vp2,[],3);

% Five-point central differences
ru = (r_um2 - 8*r_um1 + 8*r_up1 - r_up2) ./ (12*hu);
rv = (r_vm2 - 8*r_vm1 + 8*r_vp1 - r_vp2) ./ (12*hv);

% Oriented surface element
crossprod = zeros(size(ru));

crossprod(:,1) = ru(:,2).*rv(:,3) - ru(:,3).*rv(:,2);
crossprod(:,2) = ru(:,3).*rv(:,1) - ru(:,1).*rv(:,3);
crossprod(:,3) = ru(:,1).*rv(:,2) - ru(:,2).*rv(:,1);

% Evaluate vector field
Fval = F(r);
Fval = reshape(Fval,[],3);

% F dot (ru x rv)
integrand = sum(Fval .* crossprod,2);

% Integrate over parameter space
Phi = sum(integrand .* wt);

end
