function p = cyl2c(p)
rho = p(:,1);
phi = p(:,2);
z = p(:,3);

x = rho .* cos(phi);
y = rho .* sin(phi);
p = [x, y,z];
end