function p = sph2c(s)
r = s(:,1);
theta = s(:,2);
phi = s(:,3);

x = r .* sin(theta) .* cos(phi);
y = r .* sin(theta) .* sin(phi);
z = r .* cos(theta);
p = [x, y, z];
end