function s = lineCurrent(I,curve,tspan,N,rule)

if nargin < 5 || isempty(rule)
    rule = 'midpoint';
end

[pos,w,drdt] = em.quad.line(curve,tspan,N,rule);

speed = sqrt(sum(drdt.^2,2));

if any(~isfinite(speed)) || any(speed <= 0)
    error('em.src.lineCurrent:InvalidTangent', ...
        'Curve tangent must be finite and nonzero.');
end

direction = sign(tspan(2)-tspan(1));

Idl = I .* direction .* (drdt ./ speed) .* w;

s.type = 'current';
s.pos = pos;
s.w = w;
s.q = [];
s.Idl = Idl;
s.label = '';

end
