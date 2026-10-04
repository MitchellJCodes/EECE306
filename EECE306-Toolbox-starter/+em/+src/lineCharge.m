function s = lineCharge(rhoL, curve, tspan, N, rule)

% Use midpoint if rule is not provided
if nargin < 5
    rule = 'midpoint';
end

% Get positions and dl weights from quadrature
[pos,w] = em.quad.line(curve,tspan,N,rule);

% Determine charge density at each position
if isa(rhoL,'function_handle')
    q = rhoL(pos);
else
    q = rhoL * ones(size(pos,1),1);
end

% Build source structure
s.type = 'charge';
s.pos = pos;
s.w = w;
s.q = q;

end