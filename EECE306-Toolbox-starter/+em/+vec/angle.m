function theta = angle(A,B)
dotAB = sum(A .* B, 2);
magA = em.vec.mag(A);
magB = em.vec.mag(B);
x = dotAB ./ (magA .* magB);
% containing cosine boundaries to [-1, 1]
x = max(-1, min(1,x));

theta = acos(x);
end
% Computes the dot product of corresponding rows of A and B, then uses
% acos of the normalized dot product to find the angle in radians.
% Inputs A and B are Nx3 arrays and output theta is Nx1.