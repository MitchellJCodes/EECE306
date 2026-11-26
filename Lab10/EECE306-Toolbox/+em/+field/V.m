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
%
%   Reference: V = 0 at infinity. This reference is appropriate for a
%   finite-extent source; it is not finite for ideal infinite line/surface
%   distributions. A heuristic warning is issued when the represented
%   source extent is much larger than typical observation-to-source distances.

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

% A single discretized source cannot reveal how its extent changes between
% N and 2N elements. This dimensionless diagnostic instead flags a source
% whose represented bounding-box extent is at least 50 times the median
% nearest-element distance across the requested observation points. Using the
% median prevents one grid point near a point charge from triggering a false
% alarm for an otherwise compact source. It is deliberately a warning, not
% a validity test:
% compare fixed-extent N/2N discretizations when checking convergence.
if ~isempty(r_i) && ~isempty(r)
    extent = max(max(s.pos,[],1) - min(s.pos,[],1));
    nearestByObservation = min(Rmag,[],2);
    distanceScale = median(nearestByObservation);
    if distanceScale > 0 && extent / distanceScale >= 50
        warning('em.field.V:LargeSourceExtent', ...
            ['The represented source extent is at least 50 times the median ', ...
             'nearest-element distance for these observation points. Check ', ...
             'finite extent; V = 0 at infinity is not valid for an ', ...
             'ideal unbounded source.']);
    end
end

qw = reshape(s.q .* s.w, 1, []);

V_each = (1/(4*pi*em.const.eps0())) .* qw ./ Rmag;

V_sum = sum(V_each, 2);

V_r = reshape(V_sum, N, 1);

end
