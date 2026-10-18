function s = c2sph(r)
%C2SPH Cartesian to Spherical Coordinate Point Conversion
%   S = EM.COORD.C2SPH(R) returns the corresponding Nx3 spherical 
%   coordinate point s of an Nx3 Cartesian coordinate point r.
%
%   Example:
%   em.coord.c2sph([-1, 1, sqrt(6)]) % returns [2.8284, 0.5236, 2.3562]

% Get coordinates
x = r(:,1);
y = r(:,2);
z = r(:,3);

% Calculate spherical components
radial = sqrt(x.^2 + y.^2 + z.^2);
theta = acos(z ./radial);
for k = 1 : size(radial)
    if isnan(theta(k)) || isinf(theta(k))
        theta(k) = 0;
    end
end
phi = mod(atan2(y,x), 2*pi);

s = [radial, theta, phi];
end
