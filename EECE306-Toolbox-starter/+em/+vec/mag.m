function m = mag(A)
m = sqrt(sum(A.^2, 2));
end
% Intended for Nx3 array vecotrs to return the magnitude of each vector as
% Nx1