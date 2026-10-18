function d = div(F, r, h)
%DIV Divergence Operator
%   D = EM.OP.DIV(F, R, H) returns the Nx3 divergence of a vector field.
%
%   Example: s  = em.src.pointCharge(Q, [0 0 0]);
%   h = 0.05;
%   FD = @(r) em.field.D(s, r);
%   em.op.div(FD, [6 2 8], h); % returns 6.3128e-18

% Default tolerance for h
if nargin < 3
    h = 1e-5;
end

% Posistion vector safeguard
if ~isequal(size(r), [1, 3])
    error('ERROR: r must be a 1x3 position')
end

% Positive scalar constraint for h
if ~isscalar(h) || h <= 0
    error('ERROR: h must be a positive scalar')
end

% Initialize d for the loop
d = 0;

% Unit Matrix
e = diag([1 1 1]);

% Central finite difference and summation
for k = 1:3
    Fplus = F(r + h*e(k, :));
    Fminus = F(r - h*e(k, :));
    d = d + (Fplus(:, k) - Fminus(:, k))/(2*h);
end  

end

