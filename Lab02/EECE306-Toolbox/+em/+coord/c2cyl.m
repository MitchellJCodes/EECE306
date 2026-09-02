function p = c2cyl(r)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x = r(1);
y = r(2);
z = r(3);
rho = sqrt(x .^ 2 + y.^ 2);
phi = atan2(y, x);
phi = mod(phi, 2*pi);
p = [rho, phi, z];
end