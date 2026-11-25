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

% Central-difference step sizes
hu = 1e-6*(bu-au);
hv = 1e-6*(bv-av);

% Number of surface observation points
M = numel(U);

% Preallocate surface positions
pos = zeros(M,3);

pos_u_plus  = zeros(M,3);
pos_u_minus = zeros(M,3);

pos_v_plus  = zeros(M,3);
pos_v_minus = zeros(M,3);

% Evaluate surface using scalar inputs
for k = 1:M

    % Cartesian positions on surface
    pos(k,:) = surf(U(k),V(k));

    % Perturbations with respect to u
    pos_u_plus(k,:)  = surf(U(k)+hu,V(k));
    pos_u_minus(k,:) = surf(U(k)-hu,V(k));

    % Perturbations with respect to v
    pos_v_plus(k,:)  = surf(U(k),V(k)+hv);
    pos_v_minus(k,:) = surf(U(k),V(k)-hv);

end

% Partial derivative with respect to u
drdu = (pos_u_plus-pos_u_minus)./(2*hu);

% Partial derivative with respect to v
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
