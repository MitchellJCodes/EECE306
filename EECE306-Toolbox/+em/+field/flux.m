function Phi = flux(F, surf, uspan, vspan, Nu, Nv, rule)
%FLUX Flux Integral
%   PHI = EM.FIELD.FLUX(F, SURF, USPAN, VSPAN, NU, NV, RULE) returns the
%   scalar flux passing though a parametrized surface.
%
%   See also EM.QUAD.NODES.
%
%   Example: s = em.src.pointCharge(2e-4, [0 0 0]);
%   F = @(r) em.field.D(s, r);
%   surf = @(th,ph) 0.5.*[sin(th).*cos(ph) sin(th).*sin(ph) cos(th)];
%   em.field.flux(@(r) em.field.D(s, r), surf, [0 pi], [0 2*pi], 60, 120);
%   % returns 2.0002e-04
%   
%   NOTE: A backwards orientation will return a negative flux. Ex:
%   s = em.src.pointCharge(2e-4, [0 0 0]);
%   F = @(r) em.field.D(s, r);
%   surf = @(ph,th) 0.5.*[sin(th).*cos(ph) sin(th).*sin(ph) cos(th)];
%   em.field.flux(@(r) em.field.D(s, r), surf, [0 2*pi], [0 pi], 120, 60);
%   % returns -2.0002e-04

% Default rule is midpoint rule
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

