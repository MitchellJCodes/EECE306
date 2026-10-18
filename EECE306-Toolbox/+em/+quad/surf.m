function [pos, w] = surf(surface, uspan, vspan, Nu, Nv, rule)
%SURF surf-Integral Geometry Wrapper
%   [POS, W] = EM.QUD.SURF(SURFACE, USPAN, VSPAN, NU, NV, RULE) returns Mx3
%   element positions and Mx1 weights with units of area (M = NU*NV).
%
%   See also NODES.
%
%   Example: sph = @(th,ph)[sin(th).*cos(ph) sin(th).*sin(ph) cos(th)];
%   [pos, w] = em.quad.surf(sph, [0 pi], [0 2*pi], 2, 4); % returns two
%   vector quantities:
%   pos = [ 0.5000    0.5000    0.7071    w = [1.7447
%          -0.5000    0.5000    0.7071         1.7447
%          -0.5000   -0.5000    0.7071         1.7447
%           0.5000   -0.5000    0.7071         1.7447
%           0.5000    0.5000   -0.7071         1.7447
%          -0.5000    0.5000   -0.7071         1.7447
%          -0.5000   -0.5000   -0.7071         1.7447
%           0.5000   -0.5000   -0.7071]        1.7447]
%   
%   See also NODES, EM.VEC.MAG.


% Rule defaults to midpoint rule if omitted
if nargin < 6 || isempty(rule)
    rule = 'midpoint';
end

% Span of u
au = uspan(1);
bu = uspan(2);

% Span of v
av = vspan(1);
bv = vspan(2);

% Nodes and weights of each parameter
[u, wu] = em.quad.nodes(au, bu, Nu, rule);
[v, wv] = em.quad.nodes(av, bv, Nv, rule);

% 2D grid of nodes, flattened into column vectors
[U, V] = ndgrid(u,v);
U = U(:);
V = V(:);

% 2D weight grid, flattened
[WU,WV] = ndgrid(wu,wv);
w_node = abs(WU(:).*WV(:));

% Positions of the surf elements
pos = surface(U, V);
pos = reshape(pos, [], 3);

% Steps/widths of intervals
hu = eps^(1/5) .* (1 + abs(U));
hv = eps^(1/5) .* (1 + abs(V));

% Five Point Centred Difference of U
r_um2 = surface(U - 2*hu, V);
r_um1 = surface(U - hu, V);
r_up1 = surface(U + hu, V);
r_up2 = surface(U + 2*hu, V);

% Five Point Centred Difference of V
r_vm2 = surface(U, V - 2*hv);
r_vm1 = surface(U, V - hv);
r_vp1 = surface(U, V + hv);
r_vp2 = surface(U, V + 2*hv);

% Reshaping the differences into Nx3 arrays
% Reshaping U
r_um2 = reshape(r_um2,[],3);
r_um1 = reshape(r_um1,[],3);
r_up1 = reshape(r_up1,[],3);
r_up2 = reshape(r_up2,[],3);
% Reshaping V
r_vm2 = reshape(r_vm2,[],3);
r_vm1 = reshape(r_vm1,[],3);
r_vp1 = reshape(r_vp1,[],3);
r_vp2 = reshape(r_vp2,[],3);

% Partial Derivativea
ru = (r_um2 - 8*r_um1 + 8*r_up1 - r_up2) ./ (12*hu);
rv = (r_vm2 - 8*r_vm1 + 8*r_vp1 - r_vp2) ./ (12*hv);

% surf Jacobian
crossprod = zeros(size(ru));
crossprod(:,1) = ru(:,2).*rv(:,3) - ru(:,3).*rv(:,2);
crossprod(:,2) = ru(:,3).*rv(:,1) - ru(:,1).*rv(:,3);
crossprod(:,3) = ru(:,1).*rv(:,2) - ru(:,2).*rv(:,1);
J = sqrt(sum(crossprod.^2,2));

% Area weights
w = J .* w_node;
w = w(:);
end

