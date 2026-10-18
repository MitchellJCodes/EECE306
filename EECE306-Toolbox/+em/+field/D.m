function Dv = D(s, r, epsr)
%D Eletric Flux Density Function
%   DV = EM.FIELD.D(S, R, EPSR) returns the Nx3 flux density produced by a charge
%   source or a set of charge sources.
%
%   See also EM.CONST.EPS0, EM.FIELD.E, and EM.SRC.POINTCHARGE.
%
%   Example: s = em.src.pointCharge(4e-2, [1 1 0]);
%   em.field.D(s, [1 2 3]) % returns [0 1.0066e-04 3.0198e-04]

% Default epsr is 1
if nargin < 3
    epsr = 1;
end

% Electric Field
E = em.field.E(s,r);

% Constant of Permittivity
eps0 = em.const.eps0();

% Flux Calculation
Dv = epsr .* eps0 .* E;

end

