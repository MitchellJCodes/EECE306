function E_r = E(s, r)
% number of observation points
N = size(r,1);

k = 1 / (4*pi*em.const.eps0());

r_i = s.pos;

% N x M X 3 displacement vectors
r = reshape(r, [], 1, 3);
r_i = reshape(r_i, 1, [], 3);
R = r - r_i;

Rmag = sqrt(sum(R.^2, 3));
if any(Rmag(:) < 1e-12)
    error('Observation point is too close to a source element');
end

qw = reshape(s.q .* s.w, 1, [], 1);
% the squeeze function reshapes the answer into Nx3 format
E_each = k .* qw .* R ./ (Rmag.^3);
E_sum = sum(E_each, 2);
E_r = reshape(E_sum, N, 3);

end
