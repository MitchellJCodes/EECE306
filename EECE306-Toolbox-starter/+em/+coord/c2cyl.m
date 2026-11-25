function p = c2cyl(r)
x = r(:,1);
y = r(:,2);
z = r(:,3); 

rho = sqrt(x.^2 + y.^2);
phi = mod(atan2(y, x), 2*pi);
p = [rho,phi,z];
end