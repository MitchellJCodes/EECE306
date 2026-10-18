function s = surfCharge(rhoS, surf, uspan, vspan, Nu, Nv, rule)
%SURFCHARGE Surface charge Struct
%   S = SURFCHARGE(RHOS, SURF, USPAN, VSPAN, NU, NV, RULE) returns a struct
%   with a discretized list of the differential elements of a surface 
%   charge.
%
%   Example: rhoS = @(pos) 5;
%   surf = @(u,v)[u .* cos(v) u .* sin(v) zeros(size(u))];
%   em.src.surfCharge(rhoS, surf, [-1, 1], [-1, 1], 20, 20) % return struct
%   with fields:
%
%   type: 'charge'
%   pos:  [400×3 double]
%   w:    [400×1 double]
%   q:    5
%
%   See also EM.QUAD.SURF.

if nargin < 7 || isempty(rule)
    rule = 'midpoint';
end

% Positions and area weights
[pos, w] = em.quad.surf(surf, uspan, vspan, Nu, Nv, rule);

% Charge density at each position
if isa(rhoS, 'function_handle')
    q = rhoS(pos);
else
    q = rhoS * ones(size(pos, 1), 1);
end

q = q(:);
w = w(:);

if numel(q) ~= Nu*Nv
    error('rhoS must return one value per quadrature point');
end

% Source structure
s.type = 'charge';
s.pos = pos;
s.w = w;
s.q = q;
s.Idl = [];
s.label = '';
end

