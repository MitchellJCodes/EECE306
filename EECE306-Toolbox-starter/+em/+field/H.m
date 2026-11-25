function Hv = H(s, r)
% This function computes the magnetic field of the 3 dimenisional source
N = size(r,1);

c = 1 / (4*pi);

r_i = s.pos;

r = reshape(r, [], 1, 3);
r_i = reshape(r_i, 1, [], 3);
R = r - r_i;

Rmag = sqrt(sum(R.^2, 3));

if any(Rmag(:) < 1e-12)
    error('Observation point is too close to a source element');
end

Idl = reshape(s.Idl, 1, [], 3);

Cx = Idl(:,:,2).*R(:,:,3) - Idl(:,:,3).*R(:,:,2);

Cy = Idl(:,:,3).*R(:,:,1) - Idl(:,:,1).*R(:,:,3);

Cz = Idl(:,:,1).*R(:,:,2) - Idl(:,:,2).*R(:,:,1);

crossProduct = cat(3,Cx,Cy,Cz);

H_each = c .* crossProduct ./ (Rmag.^3);

Hv = reshape(sum(H_each, 2), N, 3);

end
