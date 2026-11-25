
function C = circulation(F, curve, tspan, N, rule)
%CIRCULATION Calculate the line integral of a vector field.
%   C = integral of F dot dl along a parameterized curve.

% Use midpoint if rule is not provided
if nargin < 5
    rule = 'midpoint';
end

% Integration limits
a = tspan(1);
b = tspan(2);

% Get quadrature nodes and weights
[t,w] = em.quad.nodes(a,b,N,rule);

% Calculate positions along the curve
pos = curve(t);

% Approximate tangent vectors using central differences
h = 1e-6*(b-a);

pos_plus = curve(t+h);
pos_minus = curve(t-h);

drdt = (pos_plus-pos_minus)./(2*h);

% Evaluate vector field at observation points
Fval = F(pos);

% Dot product F dot dr/dt
integrand = sum(Fval.*drdt,2);

% Integrate along the curve
C = sum(integrand.*w);

end
