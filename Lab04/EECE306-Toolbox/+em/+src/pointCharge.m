function s = pointCharge(Q, r0)

if isequal(size(r0), [1 3])
    % The '.' accesses or creates a field within the struct.
    s.type = 'charge';
    s.pos = r0;
    s.w = 1;
    s.q = Q;

else
    5
    error('r0 must be a 1x3 position');
end
