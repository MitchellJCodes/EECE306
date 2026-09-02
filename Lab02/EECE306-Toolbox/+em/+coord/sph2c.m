function r = sph2c(s)
%SPH2C Summary of this function goes here
%   Detailed explanation goes here
radial = s(1);
theta = s(2);
phi = s(3);
x = radial .* sin(theta) .* cos(phi);
y = radial .* sin(theta) .* sin(phi);
z = radial .* cos(theta);
r = [x, y, z];
end

