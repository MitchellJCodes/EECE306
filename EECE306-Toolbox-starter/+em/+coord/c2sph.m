function p = c2sph(r)
x = r(:,1);
y = r(:,2);
z = r(:,3);

r_dis = sqrt(x.^2 + y.^2 + z.^2);
theta = acos(z ./r_dis);
phi = mod(atan2(y,x), 2*pi);
p = [r_dis, theta, phi];
end