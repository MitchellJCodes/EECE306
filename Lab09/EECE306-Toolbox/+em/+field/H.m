function H_r = H(s,r)
%H Magnetic field intensity from a discretized current source (A/m).
%   H(r) = (1/(4*pi))*sum_i cross(Idl_i,r-r_i)/|r-r_i|^3.
if ~isstruct(s) || ~isfield(s,'type') || ~strcmp(s.type,'current') || ...
        ~isfield(s,'pos') || ~isfield(s,'Idl')
    error('em.field.H:InvalidSource', ...
        'Input must be a current source with pos and Idl fields.');
end
if size(r,2) ~= 3
    error('Observation points r must be N-by-3.');
end
N = size(r,1);
M = size(s.pos,1);
if size(s.pos,2) ~= 3 || ~isequal(size(s.Idl),[M 3])
    error('Current source pos and Idl must both be M-by-3.');
end
if M == 0
    H_r = zeros(N,3);
    return
end
% Vectorize over observation points and source elements.
ro = reshape(r,[],1,3);
ri = reshape(s.pos,1,[],3);
R = ro - ri;
Rmag = sqrt(sum(R.^2,3));
if any(Rmag(:) < 1e-12)
    error('Observation point is too close to a source element');
end
Idl = reshape(s.Idl,1,[],3);
C = cross(repmat(Idl,N,1),R,3);
H_each = C ./ (Rmag.^3);
H_r = reshape(sum(H_each,2)./(4*pi),N,3);
end
