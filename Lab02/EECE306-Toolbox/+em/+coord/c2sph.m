function s = c2sph(r)
%C2SPH Summary of this function goes here
%   Detailed explanation goes here
x = r(1);
y = r(2);
z = r(3);
radial = sqrt(x .^ 2 + y .^ 2 + z.^ 2);
theta = acos(z ./ radial);
phi = atan2(y, x);
s = [radial, theta, phi];
end

