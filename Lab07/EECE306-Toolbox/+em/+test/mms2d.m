function rho = mms2d(Vexact, epsr)
%MMS2D Manufacture the source density for a chosen 2-D potential.
%   rho = em.test.mms2d(Vexact,epsr) returns a function handle rho(X,Y)
%   for which Vexact approximately satisfies
%       Laplacian(Vexact) = -rho/(eps0*epsr).
%   Vexact is a function handle Vexact(X,Y). epsr may be a positive scalar
%   or a function handle epsr(X,Y). For a scalar epsr, the source is
%   -eps0*epsr*Laplacian(Vexact), estimated by centred differences.
%
%   For spatially varying epsr, this computes -div(eps0*epsr*grad(Vexact))
%   using nested centred differences. This is a numerical manufactured
%   source; use a smooth Vexact and verify grid-refinement behavior.

if ~isa(Vexact,'function_handle')
    error('Vexact must be a function handle Vexact(X,Y)');
end
if ~(isscalar(epsr) || isa(epsr,'function_handle'))
    error('epsr must be a positive scalar or a function handle epsr(X,Y)');
end
if isscalar(epsr) && (~isfinite(epsr) || epsr <= 0)
    error('epsr must be positive and finite');
end
h = eps^(1/4);
rho = @(X,Y) manufacturedSource(Vexact,epsr,X,Y,h);
end

function out = manufacturedSource(Vf,er,X,Y,h)
X = double(X); Y = double(Y);
if ~isequal(size(X),size(Y)), error('X and Y must have matching sizes'); end
if isscalar(er)
    lap = secondLaplacian(Vf,X,Y,h);
    out = -em.const.eps0()*er.*lap;
else
    % Conservative form: div(eps_r grad(V)). Fluxes are evaluated at
    % half steps, then differenced at the requested points.
    v0 = Vf(X,Y);
    qxp = er(X+h/2,Y) .* (Vf(X+h,Y)-v0)/h;
    qxm = er(X-h/2,Y) .* (v0-Vf(X-h,Y))/h;
    qyp = er(X,Y+h/2) .* (Vf(X,Y+h)-v0)/h;
    qym = er(X,Y-h/2) .* (v0-Vf(X,Y-h))/h;
    out = -em.const.eps0()*((qxp-qxm)/h + (qyp-qym)/h);
end
if any(~isfinite(out(:)))
    error('Vexact/epsr produced non-finite manufactured source values');
end
end

function lap = secondLaplacian(Vf,X,Y,h)
f0 = Vf(X,Y);
lap = (Vf(X+h,Y)-2*f0+Vf(X-h,Y))/(h^2) + ...
      (Vf(X,Y+h)-2*f0+Vf(X,Y-h))/(h^2);
end
