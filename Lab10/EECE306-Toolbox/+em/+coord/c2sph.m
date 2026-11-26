function s = c2sph(r)
%C2SPH Cartesian to Spherical Coordinate Point Conversion
%   S = EM.COORD.C2SPH(r) returns the corresponding Nx3 spherical 
%   coordinate point s of an Nx3 Cartesian coordinate point r.
%
% Example:
% em.coord.c2sph([-1, 1, sqrt(6)]) % returns [2.8284, 0.5236, 2.3562]\

x = r(:,1);
y = r(:,2);
z = r(:,3);

radial = sqrt(x.^2 + y.^2 + z.^2);
theta = zeros(size(radial));
nonzero = radial > 0;
theta(nonzero) = acos(max(-1,min(1,z(nonzero)./radial(nonzero))));
phi = mod(atan2(y,x), 2*pi);
% At the origin both angles are mathematically undefined; use zero by convention.
phi(~nonzero) = 0;
s = [radial, theta, phi];
end
