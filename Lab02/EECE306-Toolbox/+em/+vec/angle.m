function th = angle(A, B)
%ANGLE Interior Angle Function theta
% th = em.vec.angle(A, B) returns the Nx1 interior angle theta
% between two Nx3 vector inputs A and B (in radians).
%
% Example:
% em.vec.angle([0 0 1], [0, 1, 0]) % returns 1.5708 (perpendicular case)
% em.vec.angle([0 1 0], [0, 2, 0]) % returns 0 (parallel case)
% em.vec.angle([0 1 0], [0, -1, 0]) % returns 3.1416 (antiaparallel case)
%
% See also MAG.
%
% Notes: cos(theta) (cos_th) is clamped to [-1, 1] to minimize rounding 
% errors.

dot_prod = sum(A .* B, 2);
magA = em.vec.mag(A);
magB = em.vec.mag(B);
cos_th = dot_prod ./ (magA .* magB);
cos_th = max(-1, min(1, cos_th));
th = acos(cos_th);
end

