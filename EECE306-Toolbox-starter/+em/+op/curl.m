function c = curl(F,r,h)

% Default central-difference step size
if nargin < 3
    h = 1e-5;
end

% Cartesian perturbations
dx = [h 0 0];
dy = [0 h 0];
dz = [0 0 h];

% Evaluate field at perturbed positions
Fx_plus  = F(r + dx);
Fx_minus = F(r - dx);

Fy_plus  = F(r + dy);
Fy_minus = F(r - dy);

Fz_plus  = F(r + dz);
Fz_minus = F(r - dz);

% Central-difference partial derivatives
dFdx = (Fx_plus - Fx_minus)./(2*h);
dFdy = (Fy_plus - Fy_minus)./(2*h);
dFdz = (Fz_plus - Fz_minus)./(2*h);

% Cartesian curl components
Cx = dFdy(:,3) - dFdz(:,2);

Cy = dFdz(:,1) - dFdx(:,3);

Cz = dFdx(:,2) - dFdy(:,1);

% Assign output variable
c = [Cx Cy Cz];

end
