function s = volCharge(rhoV, vol, uspan, vspan, wspan, Nu, Nv, Nw, rule)
%VOLCHARGE Volume Charge Struct
%   S = VOLCHARGE(RHOS, VOL, USPAN, VSPAN, WSPAN, NU, NV, NW, RULE) 
%   returns a struct with a discretized list of the differential elements 
%   of a volume charge.
%   Example: rhoV = @(pos) 2;
%   volume = @(u, v, w) [u, v, w];
%   em.src.volCharge(rhoV, volume, [-1, 1], [-1, 1], [-1, 1], 20, 20, 20)
%   % returns a struct with fields:
%   type: 'charge'
%   pos:  [8000×3 double]
%   w:    [8000×1 double]
%   q:    2

% Default rule is midpoint rule
if nargin < 9 || isempty(rule)
    rule = 'midpoint';
end

% Position and volume weights
[pos, w] = em.quad.vol(vol, uspan, vspan, wspan, Nu, Nv, Nw, rule);

% Charge density at each point
if isa(rhoV, 'function_handle')
    q = rhoV(pos);
else
    q = rhoV * ones(size(pos, 1), 1);
end

q = q(:);
w = w(:);

if numel(q) ~= Nu*Nv*Nw
    error('rhoV must return one value per quadrature point');
end

% Source structure
s.type = 'charge';
s.pos = pos;
s.w = w;
s.q = q;
s.Idl = [];
s.label = '';

end