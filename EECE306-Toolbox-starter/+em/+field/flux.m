function Phi = flux(F, surf, uspan, vspan, Nu, Nv, rule)
if nargin < 7
    rule = 'midpoint';
end

% Parameter limits
au = uspan(1);
bu = uspan(2);

av = vspan(1);
bv = vspan(2);

% Get quadrature nodes and weights
[u, wu] = em.quad.nodes(au,bu,Nu,rule);
[v, wv] = em.quad.nodes(av,bv,Nv,rule);

% Create all combinations of u and v
[U,V] = meshgrid(u,v);

U = U(:);
V = V(:);

% Cartesian positions on surface
pos = surf(U,V);

% Central-difference step sizes
hu = 1e-6*(bu-au);
hv = 1e-6*(bv-av);

% Partial derivative with respect to u
pos_u_plus  = surf(U+hu,V);
pos_u_minus = surf(U-hu,V);

drdu = (pos_u_plus-pos_u_minus)./(2*hu);

% Partial derivative with respect to v
pos_v_plus  = surf(U,V+hv);
pos_v_minus = surf(U,V-hv);

drdv = (pos_v_plus-pos_v_minus)./(2*hv);

% Directed surface Jacobian
dS = cross(drdu,drdv,2);

% Quadrature weights
[WU,WV] = meshgrid(wu,wv);

w_node = WU(:).*WV(:);

% Evaluate vector field at surface points
Fval = F(pos);

% F dot dS
integrand = sum(Fval.*dS,2);

% Integrate over surface
Phi = sum(integrand.*w_node);

end