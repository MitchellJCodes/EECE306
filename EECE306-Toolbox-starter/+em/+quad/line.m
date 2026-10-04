function [pos, w] = line(curve, tspan, N, rule)
if nargin < 4
    rule = 'midpoint';
end

a = tspan(1);
b = tspan(2);

[t, w_node] = em.quad.nodes(a,b,N,rule);

pos = curve(t);

h = 1e-6 * (b-a);

pos_plus = curve(t+h);
pos_minus = curve(t-h);

drdt = (pos_plus - pos_minus) ./ (2*h);

drdt_mag = em.vec.mag(drdt);

w = drdt_mag .* w_node;

end

