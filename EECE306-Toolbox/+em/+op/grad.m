function g = grad(f, r, h)
%GRAD Gradient Operator
%   G = EM.FIELD.GRAD(F, R, H) returns the Nx3 gradient of a scalar field.
%
%   Example: FV = @(r) em.field.V(s, r);
%   [X, Y] = meshgrid(linspace(0.3, 2, 40), linspace(-1, 1, 40));
%   P = [X(:) Y(:) zeros(numel(X),1)];
%   E = -em.op.grad(FV, P); % returns the Nx3 eletric field computed using
%   the gradient of the potential.

% Default tolerance for h
if nargin < 3
    h = 1e-5;
end

% Posistion vector safeguard
if size(r, 2) ~= 3
    error('ERROR: r must be an N-by-3 position array')
end

% Positive scalar constraint for h
if ~isscalar(h) || h <= 0
    error('ERROR: h must be a positive scalar')
end

% Initialize g for the loop
g = zeros(size(r));

% Unit Matrix
e = diag([1 1 1]);

% Central finite difference and summation
for k = 1:3
    fplus = f(r + h*e(k, :));
    fminus = f(r - h*e(k, :));
    g(:, k) = (fplus - fminus)/(2*h);
    
end 

end

