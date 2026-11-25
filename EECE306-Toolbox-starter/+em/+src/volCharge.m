function s = volCharge(rhoV, vol, uspan, vspan, wspan, Nu, Nv, Nw, rule)

% Use midpoint if rule is not provided
if nargin < 9
    rule = 'midpoint';
end

% Get positions and volume weights from quadrature
[pos,w] = em.quad.vol(vol,uspan,vspan,wspan,Nu,Nv,Nw,rule);

% Determine charge density at each position
if isa(rhoV,'function_handle')
    q = rhoV(pos);
else
    q = rhoV * ones(size(pos,1),1);
end

% Build source structure
s.type = 'charge';
s.pos = pos;
s.w = w;
s.q = q;

end