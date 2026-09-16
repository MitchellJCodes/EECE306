function r = sph2c(s)
%SPH2C Spherical to Cartesian Coordinate Point Conversion
%   R = EM.COORD.SPH2C(s) returns the corresponding Nx3 Cartesian 
%   coordinate point r of an Nx3 spherical coordinate point s.
%
% Example:
% em.coord.sph2c([2, -5*(pi)/6, pi/6]) % returns [-0.8660, -0.5000, 1.7321]

radial = s(:,1);
theta = s(:,2);
phi = s(:,3);

x = radial .* sin(theta) .* cos(phi);
y = radial .* sin(theta) .* sin(phi);
z = radial .* cos(theta);
r = [x, y, z];
end

