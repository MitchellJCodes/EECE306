function p = c2cyl(r)
%C2CYL Cartesian to Cylindrical Coordinate Point Conversion
%   P = EM.COORD.C2CYL(r) returns the corresponding Nx3 cylindrical 
%   coordinate point p of an Nx3 Cartesian coordinate point r.
%
% Example:
% em.coord.c2cyl([1, -3, 5]) % returns [3.1623, 5.0341, 5.0000]

x = r(:,1);
y = r(:,2);
z = r(:,3); 

rho = sqrt(x.^2 + y.^2);
phi = mod(atan2(y, x), 2*pi);
p = [rho, phi, z];
end