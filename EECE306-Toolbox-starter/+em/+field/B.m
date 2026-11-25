function Bv = B(s, r, mur)
% Calculates the magnetic flux of the source 
if nargin < 3
    mur = 1;
end
H = em.field.H(s,r);

permeability = em.const.mu0();

Bv = mur.* permeability .* H;
end