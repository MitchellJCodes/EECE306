function V_r = V(s, r)
% V  Electric potential due to a charge source.
%
%   V_r = em.field.V(s,r)
%
%   s - discretized charge source with fields:
%       s.pos : source element positions
%       s.w   : source element weights
%       s.q   : charge density values
%
%   r - N-by-3 observation points
%
%   V_r - N-by-1 electric potential values

N = size(r,1);

r_i = s.pos;

% N x M x 3 displacement vectors
r = reshape(r, [], 1, 3);
r_i = reshape(r_i, 1, [], 3);
R = r - r_i;

Rmag = sqrt(sum(R.^2, 3));

if any(Rmag(:) < 1e-12)
    error('Observation point is too close to a source element');
end

qw = reshape(s.q .* s.w, 1, []);

V_each = (1/(4*pi*em.const.eps0())) .* qw ./ Rmag;

V_sum = sum(V_each, 2);

V_r = reshape(V_sum, N, 1);

end
