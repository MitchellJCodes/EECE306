function s = lineCharge(rhoL, curve, tspan, N, rule)
%LINECHARGE Line Charge Struct
%   S = EM.SRC.LINECHARGE(RHOL, CURVE, TSPAN, N, RULE) returns a struct 
%   with a discretized list of the differential elements of a line charge.
%
%   Example: rhoL = @(pos) 3;
%   curve = @(t) [cos(t), sin(t), zeros(size(t))];
%   em.src.lineCharge(rhoL, curve, [0 2*pi], 200) % returns a struct:
%   type: 'charge'
%   pos:  [200×3 double]
%   w:    [200×1 double]
%   q:    3
%
%   See also EM.QUAD.LINE.

% Default to midpoint rule
if nargin < 5 || isempty(rule)
    rule = 'midpoint';
end

% Positions and arc length weights
[pos, w] = em.quad.line(curve, tspan, N, rule);

% Charge Density at each position
if isa(rhoL, 'function_handle')
    q = rhoL(pos);
else
    q = rhoL .* ones(size(pos, 1), 1);
end

q = q(:);
w = w(:);

if numel(q) ~= N
    error('rhoL must return one value per quadrature point');
end

% Source Structure
s.type = 'charge';
s.pos = pos;
s.w = w;
s.q = q;
s.Idl = [];
s.label = '';
end

