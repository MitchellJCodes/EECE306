function Dv = D(s, r, epsr)
if nargin < 3
    epsr = 1;
end
E = em.field.E(s,r);

epsilon0 = em.const.eps0;

Dv = epsr.* epsilon0 .* E;
end