function r = cyl2c(p)
%CYL2C Cylindrical to Cartesian Coordinate Point Conversion
%   R = EM.COORD.CYL2C(p) returns the corresponding Nx3 Cartesian 
%   coordinate r of an Nx3 cylindrical coordinate p.
%
% Example:
% em.coord.cyl2c([5, pi/6, 4]) % returns [4.3301, 2.5000, 4.0000]

rho = p(:,1);
phi = p(:,2);
z = p(:,3);

x = rho .* cos(phi);
y = rho .* sin(phi);
r = [x, y, z];
end

