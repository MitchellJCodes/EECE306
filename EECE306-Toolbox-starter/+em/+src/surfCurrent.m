function s = surfCurrent(K, surf, uspan, vspan, Nu, Nv, rule)
% Use midpoint if rule is not provided
if nargin < 7 || isempty(rule)
    rule = 'midpoint';
end

% Get positions and surface-area weights from quadrature
[pos,w] = em.quad.surf(surf,uspan,vspan,Nu,Nv,rule);

% Determine charge density at each position
if isa(K,'function_handle')
    current = K(pos);
else
    current = repmat(K, size(pos,1), 1);
end

Idl = current .* w;

% Build source structure
s.type = 'current';
s.pos = pos;
s.w = w;
s.q = [];
s.Idl = Idl;
s.label = '';

end