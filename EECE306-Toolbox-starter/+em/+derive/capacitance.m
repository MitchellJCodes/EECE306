function C = capacitance(G, epsr, mask, Vfix, opts)
% Use the Poisson solver to determine the electric potential distribution
% based on the material permittivity and specified boundary conditions.
% The resulting potential is used to calculate the stored electric energy W,
% which is then used to derive capacitance using C = 2*W/V0^2.
rho = 0;
[V, ~] = em.solve.poisson2d(G, rho, epsr, mask, Vfix, opts);
W = em.derive.energyE(G, V, epsr);

voltages = Vfix(mask);
V0 = max(voltages) - min(voltages);

C = 2 * W / (V0^2);
end