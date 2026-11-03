function D_r = D(s, r, epsr)
% D  Electric flux density due to a charge source.
%
%   D_r = em.field.D(s,r,epsr)
%   epsr - relative permittivity (optional, defaults to 1)
%
%   s - discretized charge source with fields:
%       s.pos : source element positions
%       s.w   : source element weights
%       s.q   : charge density values
%
%   r - N-by-3 observation points
%
%   D_r - N-by-3 electric flux density values

if nargin < 3 || isempty(epsr)
    epsr = 1;
end
if ~isscalar(epsr) || ~isfinite(epsr) || epsr <= 0
    error('epsr must be a positive finite scalar');
end

N = size(r,1);
if size(r,2) ~= 3
    error('r must be an N-by-3 array');
end

r_i = s.pos;

% N x M x 3 displacement vectors
r = reshape(r, [], 1, 3);
r_i = reshape(r_i, 1, [], 3);
R = r - r_i;

Rmag = sqrt(sum(R.^2, 3));

if any(Rmag(:) < 1e-12)
    error('Observation point is too close to a source element');
end

qw = reshape(s.q .* s.w, 1, [], 1);

D_each = (epsr/(4*pi)) .* qw .* R ./ (Rmag.^3);

D_sum = sum(D_each, 2);

D_r = reshape(D_sum, N, 3);

end
