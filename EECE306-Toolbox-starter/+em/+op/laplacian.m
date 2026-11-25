function L = laplacian(f, r, h)
%LAPLACIAN Calculate the Laplacian of a scalar function.
% L = em.op.laplacian(f,r,h)
% Uses second-order central differences.

if nargin < 3
    h = 1e-5;
end

% Evaluate field at original points
F = f(r);

% Create shifted points for x derivative
rx_plus = r;
rx_minus = r;

rx_plus(:,1) = rx_plus(:,1) + h;
rx_minus(:,1) = rx_minus(:,1) - h;

% Evaluate field
Fx_plus = f(rx_plus);
Fx_minus = f(rx_minus);

% Second derivative with respect to x
d2Fdx2 = (Fx_plus - 2*F + Fx_minus) ./ (h^2);


% Create shifted points for y derivative
ry_plus = r;
ry_minus = r;

ry_plus(:,2) = ry_plus(:,2) + h;
ry_minus(:,2) = ry_minus(:,2) - h;

% Evaluate field
Fy_plus = f(ry_plus);
Fy_minus = f(ry_minus);

% Second derivative with respect to y
d2Fdy2 = (Fy_plus - 2*F + Fy_minus) ./ (h^2);


% Create shifted points for z derivative
rz_plus = r;
rz_minus = r;

rz_plus(:,3) = rz_plus(:,3) + h;
rz_minus(:,3) = rz_minus(:,3) - h;

% Evaluate field
Fz_plus = f(rz_plus);
Fz_minus = f(rz_minus);

% Second derivative with respect to z
d2Fdz2 = (Fz_plus - 2*F + Fz_minus) ./ (h^2);


% Laplacian
L = d2Fdx2 + d2Fdy2 + d2Fdz2;

end
