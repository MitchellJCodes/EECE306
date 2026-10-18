function s = pointCharge(Q, r0)
%POINTCHARGE Point Source Struct
%   S = EM.SRC.POINTCHARGE(Q, r0) returns a struct with a discretized list
%   of differential elements of a point charge. It contains the type,
%   position, quadrature weights, and charge density of a point source.
%
%   Example:
%   em.src.pointCharge(1e-12, [0 0 0]) returns struct with fields:
%   type: 'charge'
%   pos: [0 0 0]
%   w: 1
%   q: 1.0000e-12

if isequal(size(r0), [1 3])
    % The '.' accesses or creates a field within the struct.
    s.type = 'charge';
    s.pos = r0;
    s.w = 1;
    s.q = Q;
    s.Idl = [];
    s.label = '';

else
    error('ERROR: r0 must be a 1x3 position');
end

