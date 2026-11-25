function rho = mms2d(Vexact, epsr)
%MMS2D Manufactured source for Poisson's equation.
% rho = em.test.mms2d(Vexact, epsr)
%
% Vexact is a function handle @(X,Y)
% epsr is the uniform relative permittivity
% rho is returned as a function handle @(X,Y)

e0 = em.const.eps0();

rho = @(X,Y) source(X,Y,Vexact,epsr,e0);

end

function R = source(X,Y,Vexact,epsr,e0)

% Convert 2D coordinates to N-by-3 observation points
r = [X(:), Y(:), zeros(numel(X),1)];

% Wrap the 2D potential for the 3D Laplacian helper
f = @(p) Vexact(p(:,1),p(:,2));

% Calculate Laplacian
L = em.op.laplacian(f,r,1e-4);

% Poisson equation
R = reshape(-e0 .* epsr .* L,size(X));

end
