function C = circulation(F, curve, tspan, N, rule)
%CIRCULATION Flux Integral on Closed Path
%   C = EM.FIELD.CIRCULATION(F, CURVE, TSPAN, N, RULE) % returns the
%   scalar flux of a vector field around a curve.
%
%   See also EM.QUAD.NODES.
%
%   Example: FE = @(r) em.field.E(s, r);
%   loops = { @(t)[2.*cos(t) 2.*sin(t) zeros(size(t))], ...
%         @(t)[1+0.5.*cos(t) 0.5*sin(3*t) 0.3.*sin(t)], ...
%         @(t)[0.5.*cos(t) 1+0.8.*sin(t) 0.2.*cos(2.*t)] };
%   for k = 1:3
%       C = em.field.circulation(FE, loops{k}, [0 2*pi], 400);
%   end  % returns -4.6745e11

% Default rule is midpoint rule
if nargin < 7
    rule = 'midpoint';
end

% Intervals
ta = tspan(1);
tb = tspan(2);

% Parameter-space quadrature
[t, wt] = em.quad.nodes(ta, tb, N, rule);

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

