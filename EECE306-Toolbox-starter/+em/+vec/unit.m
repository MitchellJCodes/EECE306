function u = unit(A)
m = em.vec.mag(A);
u = A ./ m;
end
% Computes the unit vector of each row of A.
% Input A is Nx3 and output is Nx3.