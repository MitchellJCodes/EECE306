function [pos, w] = vol(volfun, uspan, vspan, wspan, Nu, Nv, Nw, rule)
if nargin < 8
    rule = 'midpoint';
end

au = uspan(1);
bu = uspan(2);

av = vspan(1);
bv = vspan(2);

aw = wspan(1);
bw = wspan(2);

[u, wu] = em.quad.nodes(au,bu,Nu,rule);
[v, wv] = em.quad.nodes(av,bv,Nv,rule);
[w, ww] = em.quad.nodes(aw,bw,Nw,rule);

[U,V, W] = meshgrid(u,v, w);
U = U(:);
V = V(:);
W = W(:);

pos = volfun(U,V,W);

hu = 1e-6 * (bu-au);

pos_u_plus = volfun(U+hu, V, W);
pos_u_minus = volfun(U-hu,V, W);

drdu = (pos_u_plus - pos_u_minus) ./ (2*hu);

hv = 1e-6 * (bv-av);

pos_v_plus = volfun(U,V+hv, W);
pos_v_minus = volfun(U,V-hv, W);

drdv = (pos_v_plus - pos_v_minus) ./ (2*hv);

hw = 1e-6 * (bw-aw);

pos_w_plus = volfun(U,V, W+hw);
pos_w_minus = volfun(U,V, W-hw);

drdw = (pos_w_plus - pos_w_minus) ./ (2*hw);

J = abs(dot(drdu,cross(drdv,drdw,2),2));

[WU, WV, WW] = meshgrid(wu, wv, ww);

w_node = WU(:) .* WV(:) .* WW(:);
w = J .* w_node;

end