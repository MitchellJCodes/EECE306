
function V = V(s,r)
%V Calculate electric potential from a charge distribution.
%   V = em.field.V(s,r)
%   Reference: V = 0 at infinity for finite sources.

% Warn when source extent grows with element count
if isfield(s,'growingExtent') && s.growingExtent
    warning('em:field:InfiniteSource', ...
        ['Source extent grows with element count. ', ...
         'Potential referenced to infinity may diverge.']);
end

% Coulomb constant
k = 1/(4*pi*em.const.eps0());

r_i = s.pos;

% N x M x 3 displacement vectors
r = reshape(r,[],1,3);
r_i = reshape(r_i,1,[],3);

R = r-r_i;

% Distance from each source element
Rmag = sqrt(sum(R.^2,3));

% Check for singularities
if any(Rmag(:) < 1e-12)
    error('Observation point is too close to a source element');
end

% Charge contribution from each source element
qw = (s.q.*s.w).';

% Potential contribution from each element
V_each = k.*qw./Rmag;

% Sum contributions to obtain N x 1 potential
V = sum(V_each,2);

end
