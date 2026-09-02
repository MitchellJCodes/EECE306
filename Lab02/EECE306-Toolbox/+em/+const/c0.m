function c = c0()
%C0 Speed Of Light Constant c
% c = em.const.c0() returns the scalar value of the speed of light.
%
% Example:
% em.const.c0() % returns 2.9979e+08
%
% See also EPS0, MU0.
%
% Notes: This constant is computed from em.const.eps0() and em.const.mu0().
eps0 = em.const.eps0();
mu0 = em.const.mu0();
c = 1/sqrt(eps0 * mu0);
end