function g = grad(f, r, h)
if nargin < 3
    h = 1e-5;
end

% Create shifted points for x derivative
rx_plus = r;
rx_minus = r;

rx_plus(:,1) = rx_plus(:,1) + h;
rx_minus(:,1) = rx_minus(:,1) - h;

% Evaluate field
Fx_plus = f(rx_plus);
Fx_minus = f(rx_minus);

% dFx/dx
dFxdx = (Fx_plus - Fx_minus) ./ (2*h);


% Create shifted points for y derivative
ry_plus = r;
ry_minus = r;

ry_plus(:,2) = ry_plus(:,2) + h;
ry_minus(:,2) = ry_minus(:,2) - h;

% Evaluate field
Fy_plus = f(ry_plus);
Fy_minus = f(ry_minus);

% dFy/dy
dFydy = (Fy_plus - Fy_minus) ./ (2*h);


% Create shifted points for z derivative
rz_plus = r;
rz_minus = r;

rz_plus(:,3) = rz_plus(:,3) + h;
rz_minus(:,3) = rz_minus(:,3) - h;

% Evaluate field
Fz_plus = f(rz_plus);
Fz_minus = f(rz_minus);

% dFz/dz
dFzdz = (Fz_plus - Fz_minus) ./ (2*h);


% Gradient
g = [dFxdx, dFydy, dFzdz];

end