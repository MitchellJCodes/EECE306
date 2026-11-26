function C = curl(F, r0, h)
%CURL Compute the numerical curl of a vector field using central differences.
%
%   C = em.op.curl(F,r0,h)
%
%   F   - vector field function handle accepting N-by-3 points and
%         returning N-by-3 vectors
%   r0  - 1-by-3 or N-by-3 observation point(s)
%   h   - positive central-difference step size (default 1e-5)
%
%   C   - 1-by-3 or N-by-3 curl of F at r0
%
%   The Cartesian components are
%     curl(F)_x = dFz/dy - dFy/dz
%     curl(F)_y = dFx/dz - dFz/dx
%     curl(F)_z = dFy/dx - dFx/dy

if nargin < 3 || isempty(h)
    h = 1e-5;
end

if ~isa(F, 'function_handle')
    error('em.op.curl:InvalidField', 'F must be a function handle.');
end
if ~isnumeric(r0) || ~isreal(r0) || size(r0,2) ~= 3 || ...
        isempty(r0) || any(~isfinite(r0(:)))
    error('em.op.curl:InvalidPoints', ...
        'r0 must be a nonempty finite N-by-3 array.');
end
if ~isscalar(h) || ~isreal(h) || ~isfinite(h) || h <= 0
    error('em.op.curl:InvalidStep', 'h must be a positive finite scalar.');
end

% Evaluate the field at the six displaced point sets.
FxP = evalVectorField(F, r0 + [h 0 0]);
FxM = evalVectorField(F, r0 - [h 0 0]);
FyP = evalVectorField(F, r0 + [0 h 0]);
FyM = evalVectorField(F, r0 - [0 h 0]);
FzP = evalVectorField(F, r0 + [0 0 h]);
FzM = evalVectorField(F, r0 - [0 0 h]);

% Central differences, with the cyclic Cartesian curl components.
dFz_dy = (FyP(:,3) - FyM(:,3)) ./ (2*h);
dFy_dz = (FzP(:,2) - FzM(:,2)) ./ (2*h);

dFx_dz = (FzP(:,1) - FzM(:,1)) ./ (2*h);
dFz_dx = (FxP(:,3) - FxM(:,3)) ./ (2*h);

dFy_dx = (FxP(:,2) - FxM(:,2)) ./ (2*h);
dFx_dy = (FyP(:,1) - FyM(:,1)) ./ (2*h);

C = [dFz_dy - dFy_dz, ...
     dFx_dz - dFz_dx, ...
     dFy_dx - dFx_dy];
end

function value = evalVectorField(F, points)
value = F(points);
n = size(points,1);
if ~isnumeric(value) || ~isreal(value) || numel(value) ~= 3*n
    error('em.op.curl:InvalidFieldOutput', ...
        'F must return one real three-component vector per point.');
end
value = reshape(value, n, 3);
if any(~isfinite(value(:)))
    error('em.op.curl:NonfiniteField', ...
        'F returned a nonfinite value at a displaced observation point.');
end
end
