function s = surfCharge(rhoS, surf, uspan, vspan, Nu, Nv, rule)

% Use midpoint if rule is not provided
if nargin < 7
    rule = 'midpoint';
end

% Get positions and surface-area weights from quadrature
[pos,w] = em.quad.surf(surf,uspan,vspan,Nu,Nv,rule);

% Determine charge density at each position
if isa(rhoS,'function_handle')
    q = rhoS(pos);
else
    q = rhoS * ones(size(pos,1),1);
end

% Build source structure
s.type = 'charge';
s.pos = pos;
s.w = w;
s.q = q;

end