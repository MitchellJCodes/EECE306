function g = grad(f, r0, h)
%GRAD Compute the numerical gradient of a scalar field.
%
%   g = em.op.grad(f,r0,h)
%
%   f   - scalar field function handle
%   r0  - 1-by-3 observation point
%   h   - finite-difference step size
%
%   g   - 1-by-3 gradient of f at r0

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

% Evaluate the scalar field
fx_plus  = f(rxp);
fx_minus = f(rxm);

fy_plus  = f(ryp);
fy_minus = f(rym);

fz_plus  = f(rzp);
fz_minus = f(rzm);

% Central finite differences
df_dx = (fx_plus - fx_minus) / (2*h);
df_dy = (fy_plus - fy_minus) / (2*h);
df_dz = (fz_plus - fz_minus) / (2*h);

% Gradient
g = [df_dx df_dy df_dz];

end
