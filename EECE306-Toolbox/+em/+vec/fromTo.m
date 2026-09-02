function d = fromTo(P, Q)
%FROMTO Displacement Vector Function
%  d = em.vec.fromTo(P, Q) returns the Nx3 displacement vector d (= Q - P) 
% between two Nx3 vectors P and Q
%
% Example:
% em.vec.fromTo([1 2 3], [4, 5, 6]) % returns [3 3 3]

d = Q - P;
end