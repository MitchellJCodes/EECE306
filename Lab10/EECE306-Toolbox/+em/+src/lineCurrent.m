function s = lineCurrent(I,curve,tspan,N,rule)
%LINECURRENT Discretize a filamentary current along a parametrized curve
%   CURVE must return positions in R^3
if nargin < 5 || isempty(rule), rule = 'midpoint'; end
[pos,w,drdt] = em.quad.line(curve,tspan,N,rule);
speed = sqrt(sum(drdt.^2,2));
if any(~isfinite(speed)) || any(speed <= 0)
    error('em.src.lineCurrent:InvalidTangent', ...
        'Curve tangent must be finite and nonzero at every quadrature node.');
end
% w is the positive physical line-element weight
direction = sign(tspan(2)-tspan(1));
Idl = I .* direction .* (drdt ./ speed) .* w;
s.type = 'current';
s.pos = pos;
s.w = w;
s.q = [];
s.Idl = Idl;
s.label = '';
end
