function eta = eta0()
%ETA0 Intrinsic Impedance Of Free Space eta (eta0)
% eta = em.const.eta0() returns the scalar value of the impedance of free
% space.
%
% Example:
% em.const.eta0() % returns 376.7303
%
% See also EPS0, MU0.
%
% Notes: This constant is computed from em.const.eps0() and em.const.mu0().
eps0 = em.const.eps0();
mu0 = em.const.mu0();
eta = sqrt(mu0 / eps0);
end