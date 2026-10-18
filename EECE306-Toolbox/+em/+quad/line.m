function [pos, w] = line(curve, tspan, N, rule)
%LINE Line-Integral Geometry Wrapper
%   [POS, W] = EM.QUAD.LINE(CURVE, TSPAN, N, RULE) returns the Nx3 element
%   positions POS and Nx1 weights W with units of length.
%
%   Example: circ = @(t)[cos(t) sin(t) zeros(size(t))];
%   [pos, w] = em.quad.line(circ, [0 2*pi], 6) % returns two vectors:
%   pos = [ 0.8660        0.5000   0.0000        w = [ 1.0472
%           6.1232e-17    1.0000   0.0000              1.0472
%          -0.8660        0.5000   0.0000              1.0472
%          -0.8660       -0.5000   0.0000              1.0472
%          -1.8370e-16   -1.5000   0.0000              1.0472
%           0.8660       -0.5000   0.0000 ]            1.0472]
%
%   See also EM.VEC.MAG, NODES.
%   
%   NOTE: curve is a function handle @(t) -> 1x3 or Nx3.

if nargin < 4 || isempty(rule)
    rule = 'midpoint';
end

% tspan (= [a, b]) is the parameter domain of the curve
a = tspan(1); % Left bound
b = tspan(2); % Right Bound

[t, w_node] = em.quad.nodes(a, b, N, rule); % Node and weight vectors
pos = curve(t); % Coordinates of each line element
pos = reshape(pos,[],3); % Positions must be Nx3

% Step
h = eps^(1/5) .* (1 + abs(t));

% Five-point centered difference
rm2 = curve(t - 2*h);
rm1 = curve(t - h);
rp1 = curve(t + h);
rp2 = curve(t + 2*h);

% Reshape each difference into Nx3 arrays
rm2 = reshape(rm2,[],3);
rm1 = reshape(rm1,[],3);
rp1 = reshape(rp1,[],3);
rp2 = reshape(rp2,[],3);

% Physical line element: ds = |dr/dt| dt, where |dr/dt| = speed

% Derivative and speed
drdt = (rm2 - 8*rm1 + 8*rp1 - rp2) ./ (12*h);
speed = sqrt(sum(drdt.^2,2));

% Arc length weights
% abs(wt) ensures the length is always positive
w = abs(w_node).*speed;
w = w(:);

end







