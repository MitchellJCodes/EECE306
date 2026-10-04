function d = div(F, r, h)
if nargin < 3
    h = 1e-5;
end

% Create shifted points for x derivative
rx_plus = r;
rx_minus = r;

rx_plus(:,1) = rx_plus(:,1) + h;
rx_minus(:,1) = rx_minus(:,1) - h;

% Evaluate field
Fx_plus = F(rx_plus);
Fx_minus = F(rx_minus);

% dFx/dx
dFxdx = (Fx_plus(:,1) - Fx_minus(:,1)) ./ (2*h);


% Create shifted points for y derivative
ry_plus = r;
ry_minus = r;

ry_plus(:,2) = ry_plus(:,2) + h;
ry_minus(:,2) = ry_minus(:,2) - h;

% Evaluate field
Fy_plus = F(ry_plus);
Fy_minus = F(ry_minus);

% dFy/dy
dFydy = (Fy_plus(:,2) - Fy_minus(:,2)) ./ (2*h);


% Create shifted points for z derivative
rz_plus = r;
rz_minus = r;

rz_plus(:,3) = rz_plus(:,3) + h;
rz_minus(:,3) = rz_minus(:,3) - h;

% Evaluate field
Fz_plus = F(rz_plus);
Fz_minus = F(rz_minus);

% dFz/dz
dFzdz = (Fz_plus(:,3) - Fz_minus(:,3)) ./ (2*h);


% Divergence
d = dFxdx + dFydy + dFzdz;

end