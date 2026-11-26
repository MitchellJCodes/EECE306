function d = div(F, r0, h)
%DIV Compute the numerical divergence of a vector field.
%
%   d = em.op.div(F,r0,h)
%
%   F  - vector field function handle
%   r0 - 1-by-3 observation point
%   h  - finite-difference step size
%
%   d  - divergence of F at r0

if nargin < 3 || isempty(h)
    h = 1e-5;
end

if ~isequal(size(r0), [1 3])
    error('r0 must be a 1x3 position');
end

if ~isscalar(h) || h <= 0
    error('h must be a positive scalar');
end

% Points used for the central finite differences
rxp = r0 + [h 0 0];
rxm = r0 - [h 0 0];

ryp = r0 + [0 h 0];
rym = r0 - [0 h 0];

rzp = r0 + [0 0 h];
rzm = r0 - [0 0 h];

% Evaluate the vector field
Fx_plus  = F(rxp);
Fx_minus = F(rxm);

Fy_plus  = F(ryp);
Fy_minus = F(rym);

Fz_plus  = F(rzp);
Fz_minus = F(rzm);

% Central finite differences
dFx_dx = (Fx_plus(1) - Fx_minus(1)) / (2*h);
dFy_dy = (Fy_plus(2) - Fy_minus(2)) / (2*h);
dFz_dz = (Fz_plus(3) - Fz_minus(3)) / (2*h);

% Divergence
d = dFx_dx + dFy_dy + dFz_dz;

end
