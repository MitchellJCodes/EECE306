function [pos, w] = surf(surfFun, uspan, vspan, Nu, Nv, rule)
if nargin < 6
    rule = 'midpoint';
end

au = uspan(1);
bu = uspan(2);

av = vspan(1);
bv = vspan(2);

[u, wu] = em.quad.nodes(au,bu,Nu,rule);
[v, wv] = em.quad.nodes(av,bv,Nv,rule);

[U,V] = meshgrid(u,v);
U = U(:);
V = V(:);

pos = surfFun(U,V);

hu = 1e-6 * (bu-au);

pos_u_plus = surfFun(U+hu, V);
pos_u_minus = surfFun(U-hu,V);

drdu = (pos_u_plus - pos_u_minus) ./ (2*hu);

hv = 1e-6 * (bv-av);

pos_v_plus = surfFun(U,V+hv);
pos_v_minus = surfFun(U,V-hv);

drdv = (pos_v_plus - pos_v_minus) ./ (2*hv);

J = em.vec.mag(cross(drdu,drdv,2));

[WU, WV] = meshgrid(wu, wv);

w_node = WU(:) .* WV(:);
w = J .* w_node;

end