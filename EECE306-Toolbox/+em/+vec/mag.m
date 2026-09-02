function m = mag(A)
%MAG Vector Magnitude Function
% m = em.vec.mag(A) returns the Nx1 magnitude m of an Nx3 vector A.
%
% Example:
% em.vec.mag([3 4 0]) % returns 5

m = sqrt(sum(A.^2, 2));
end