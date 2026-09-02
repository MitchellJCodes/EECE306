function u = unit(A)
%UNIT Normalization Function u
% u = em.vec.unit(A) returns the normalized Nx3 vector u of a non-zero Nx3 
% vector A
%
% Example:
% em.vec.unit([6, 8, 0]) % returns [0.6000, 0.8000, 0]
%
% See also MAG.
%
% Notes: A zero row returns [NaN, NaN, NaN] due to division by 0.

u = A ./ em.vec.mag(A);
end

