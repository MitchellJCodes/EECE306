function V = V(s, r)
%V Electric Potential Function
%   V = EM.FIELD.V(S, R) returns the Nx1 potentials of points inside an
%   electric field.
%
%   See also EM.CONST.EPS0.
%
%   Example: Q = 2e-9;
%   radii = (1:10)';
%   s  = em.src.pointCharge(Q, [0 0 0]);
%   V = em.field.V(s, [radii zeros(10,2)]); % returns the 10x1 vector:
%   V = [17.9751
%         8.9876
%         5.9917
%         4.4938
%         3.5950
%         2.9959
%         2.5678
%         2.2469
%         1.9972
%         1.7975]
%
% NOTE: Ensure that the extent of the potential does not scale with N. If
% one of the bounds of the span is set to N (ex: span = [-N 3], N = 3), the
% potential diverges with increasing N.

% Number of sample points
N = size(r,1);

% Coulomb constant
k = 1 / (4*pi*em.const.eps0());

% Position vector
r_i = s.pos;

% N x M X 3 displacement vectors
r = reshape(r, [], 1, 3);
r_i = reshape(r_i, 1, [], 3);
R = r - r_i;

% Magnitude of displacement vector
Rmag = sqrt(sum(R.^2, 3));
tol = 1e-9;

% Singularity safeguard
[~, off_elem] = find(Rmag < tol);
if ~isempty(off_elem)
    error('ERROR: Singularity found. Offending element: #%d', off_elem(1));
end

% Growing extent test
if ~isempty(r_i) && ~isempty(r)
    extent = max(max(s.pos,[],1) - min(s.pos,[],1));
    nearestByObservation = min(Rmag,[],2);
    distanceScale = median(nearestByObservation);
    if distanceScale > 0 && extent / distanceScale >= 50
        warning(['The represented source extent is at least 50 times ' ...
            'the median nearest-element distance for these ' ...
            'observation points. Check finite extent; V = 0 at infinity' ...
            ' is not valid for an ideal unbounded source.']);
    end
end

% Total potential contributed by all points
qw = reshape(s.q .* s.w, 1, []);
V_elem = k .* qw ./ Rmag;
V_tot = sum(V_elem, 2);
V = reshape(V_tot, N, 1);

end

