function R = resistance(G, sigma, mask, Vfix, opts)
% Resistance is calculated by first solving for the electric potential
% using the steady-state conduction equation. The potential gradient
% determines the electric field, which is used with conductivity to
% calculate dissipated power. Resistance is then derived using R = V0^2/P.
rho =0;
[V, ~] = em.solve.poisson2d(G, rho, sigma, mask, Vfix, opts);

voltages = Vfix(mask);
V0 = max(voltages) - min(voltages);
% Grid Spacing
dx = G.x(2) - G.x(1);
dy = G.y(2) - G.y(1);

[dVdx, dVdy] = gradient(V, dx, dy);

E2 = dVdx.^2 + dVdy.^2;
powerDensity = sigma .* E2;
Px = trapz(G.x(:).', powerDensity, 2);
P = trapz(G.y(:), Px, 1);

R = (V0^2) / P;
end