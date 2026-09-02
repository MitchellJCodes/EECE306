function r = cyl2c(p)
%CYL2C Summary of this function goes here
%   Detailed explanation goes here
rho = p(1);
phi = p(2);
z = p(3);
x = rho .* cos(phi);
y = rho .* sin(phi);
r = [x, y, z];

end

