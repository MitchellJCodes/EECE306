function C = circulation(F, curve, interval, N)
%CIRCULATION Compute the circulation of a vector field around a curve.
%
%   C = em.field.circulation(F,curve,interval,N)
%
%   F        - vector field function handle
%   curve    - parametrized curve r(t)
%   interval - parameter interval [a b]
%   N        - number of quadrature points
%
%   C        - signed circulation of F along the curve
%
%   The circulation is
%
%       C = integral_C F . dr
%
%   evaluated using midpoint quadrature and a numerical derivative
%   of the parametrized curve.

% Parameter-space quadrature
[t,wt] = em.quad.nodes(interval(1),interval(2),N,'midpoint');

% Curve positions
r = curve(t);
r = reshape(r,[],3);

% Numerical derivative step
h = eps^(1/5) .* (1 + abs(t));

% Five-point derivative
r_m2 = curve(t - 2*h);
r_m1 = curve(t - h);
r_p1 = curve(t + h);
r_p2 = curve(t + 2*h);

% Reshape derivative results
r_m2 = reshape(r_m2,[],3);
r_m1 = reshape(r_m1,[],3);
r_p1 = reshape(r_p1,[],3);
r_p2 = reshape(r_p2,[],3);

% Five-point central difference
drdt = (r_m2 - 8*r_m1 + 8*r_p1 - r_p2) ./ (12*h);

% Evaluate vector field
Fval = F(r);
Fval = reshape(Fval,[],3);

% F dot dr/dt
integrand = sum(Fval .* drdt,2);

% Integrate over parameter space
C = sum(integrand .* wt);

end
