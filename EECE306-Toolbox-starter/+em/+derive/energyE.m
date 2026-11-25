function W = energyE(G, V, epsr)
% Defining permitivity of free space, Grid spacing, and gradient
% electric potential, while calculating electrostatic energy using the
% discrete gradient method.


eps0 = em.const.eps0();
% Grid Spacing
dx = G.x(2) - G.x(1);
dy = G.y(2) - G.y(1);

if isscalar(epsr)
    epsr = epsr * ones(size(V));
end

Ex = -(V(:,2:end) - V(:,1:end-1)) / dx;
Ey = -(V(2:end,:) - V(1:end-1,:)) / dy;

epsX = 2 .* epsr(:,1:end-1) .* epsr(:,2:end) ./ ...
       (epsr(:,1:end-1) + epsr(:,2:end));

epsY = 2 .* epsr(1:end-1,:) .* epsr(2:end,:) ./ ...
       (epsr(1:end-1,:) + epsr(2:end,:));

wx = ones(size(Ex));
wx(1,:) = 0.5;
wx(end,:) = 0.5;

wy = ones(size(Ey));
wy(:,1) = 0.5;
wy(:,end) = 0.5;

Wx = 0.5 * eps0 * sum(wx .* epsX .* Ex.^2, 'all') * dx * dy;
Wy = 0.5 * eps0 * sum(wy .* epsY .* Ey.^2, 'all') * dx * dy;

W = Wx + Wy;

end
