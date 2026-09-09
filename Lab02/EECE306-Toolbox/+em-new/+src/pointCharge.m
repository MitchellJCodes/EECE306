function s = pointCharge(Q, r0)
%POINTCHARGE Summary of this function goes here
%   Detailed explanation goes here
if isequal(size(r0), [1 3])
    s.type = 'charge';
    s.pos = r0;
    s.w = 1;
    s.q = Q;

else
    5
    error('r0 must be a 1x3 position');
end

