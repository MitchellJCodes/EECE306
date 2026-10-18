function [pos, weight] = vol(volume, uspan, vspan, wspan, Nu, Nv, Nw, rule)
%VOL Volume-Interval Geometry Wrapper
%   [POS, W] = EM.QUAD.VOL(VOL, USPAN, VSPAN, WSPAN, NU, NV, NW, RULE)
%   returns Mx3 element positions and Mx1 weights with units of volume.
%
%   See also NODES.
%
%   Example:
%   ball = @(u,v,w)[w.*sin(u).*cos(v) w.*sin(u).*sin(v) w.*cos(u)];
%   [pos, w]= em.quad.vol(ball, [0 pi], [0 2*pi], [0 1], 1, 4, 2);
%   % returns two vectors:
%   pos = [ 0.1768    0.1768    1.5308e-17    w = [ 0.1542
%          -0.1768    0.1768    1.5308e-17          0.1542
%          -0.1768   -0.1768    1.5308e-17          0.1542
%           0.1768   -0.1768    1.5308e-17          0.1542
%           0.5303    0.5303    4.5924e-17          0.1542
%          -0.5303    0.5303    4.5924e-17          0.1542
%          -0.5303   -0.5303    4.5924e-17          0.1542
%           0.5303   -0.5303    4.5924e-17 ]        0.1542 ]

% Midpoint rule is default rule
if nargin < 8
    rule = 'midpoint';
end

% Span of u
au = uspan(1);
bu = uspan(2);

% Span of v
av = vspan(1);
bv = vspan(2);

% Span of w
aw = wspan(1);
bw = wspan(2);

% Nodes and weights of each parameter
[u, wu] = em.quad.nodes(au,bu,Nu,rule);
[v, wv] = em.quad.nodes(av,bv,Nv,rule);
[w, ww] = em.quad.nodes(aw,bw,Nw,rule);

% 3D grid of nodes, flattened to column vectors
[U, V, W] = meshgrid(u, v, w);
U = U(:);
V = V(:);
W = W(:);

% Positions of the volume elements
pos = volume(U, V, W);

% 3D weight grid, flattened
[WU, WV, WW] = meshgrid(wu, wv, ww);
w_node = WU(:) .* WV(:) .* WW(:);

% Steps/widths of the intervals
hu = eps^(1/5) .* (1 + abs(U));
hv = eps^(1/5) .* (1 + abs(V));
hw = eps^(1/5) .* (1 + abs(W));

% Centered Differences with respect to U
r_um2 = volume(U-2*hu, V, W);
r_um1 = volume(U-hu, V, W);
r_up1 = volume(U+hu, V, W);
r_up2 = volume(U+2*hu, V, W);

% Centered Differences with respect to V
r_vm2 = volume(U, V-2*hv, W);
r_vm1 = volume(U, V-hv, W);
r_vp1 = volume(U, V+hv, W);
r_vp2 = volume(U, V+2*hv, W);

% Centered Differences with respect to W
r_wm2 = volume(U, V, W-2*hw);
r_wm1 = volume(U, V, W-hw);
r_wp1 = volume(U, V, W+hw);
r_wp2 = volume(U, V, W+2*hw);

% Reshape U to Nx3 array
r_um2 = reshape(r_um2,[],3);
r_um1 = reshape(r_um1,[],3);
r_up1 = reshape(r_up1,[],3);
r_up2 = reshape(r_up2,[],3);

% Reshape V to Nx3 array
r_vm2 = reshape(r_vm2,[],3);
r_vm1 = reshape(r_vm1,[],3);
r_vp1 = reshape(r_vp1,[],3);
r_vp2 = reshape(r_vp2,[],3);

% Reshape W to Nx3 array
r_wm2 = reshape(r_wm2,[],3);
r_wm1 = reshape(r_wm1,[],3);
r_wp1 = reshape(r_wp1,[],3);
r_wp2 = reshape(r_wp2,[],3);

ru = (r_um2 - 8*r_um1 + 8*r_up1 - r_up2) ./ (12*hu);
rv = (r_vm2 - 8*r_vm1 + 8*r_vp1 - r_vp2) ./ (12*hv);
rw = (r_wm2 - 8*r_wm1 + 8*r_wp1 - r_wp2) ./ (12*hw);

% Jacobian
cross_uv = zeros(size(ru));
cross_uv(:,1) = ru(:,2).*rv(:,3) - ru(:,3).*rv(:,2);
cross_uv(:,2) = ru(:,3).*rv(:,1) - ru(:,1).*rv(:,3);
cross_uv(:,3) = ru(:,1).*rv(:,2) - ru(:,2).*rv(:,1);
detJ = cross_uv(:,1).*rw(:,1) + cross_uv(:,2).*rw(:,2) + cross_uv(:,3).*rw(:,3);
J = abs(detJ);

weight = w_node .* J;
weight = weight(:);

end

