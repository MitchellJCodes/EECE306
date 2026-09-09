function quiver2(F, xlim, ylim, N, opts)

figure
% Create N evenly spaced x and y coordinates within the plot limits
x = linspace(xlim(1), xlim(2), N);
y = linspace(ylim(1), ylim(2), N);

% Create an N x N grid of x and y observation locations
[X,Y] = meshgrid(x, y);

% Convert the grid into Nx3 observation points in the z = 0 plane
r = [X(:), Y(:), zeros(numel(X),1)];

% Calculate the electric field at each observation point
E = F(r);

% Reshape the x and y electric field components to match the grid
Ex = reshape(E(:,1), size(X));
Ey = reshape(E(:,2), size(Y));

% If normalize is true, calculate the field magnitude and make
% the arrows unit length to emphasize field direction
if opts.normalize
    M = sqrt(Ex.^2 + Ey.^2);

    % Prevent division by zero
    M(M == 0) = 1;

    Ex = Ex ./ M;
    Ey = Ey ./ M;
end

quiver(X, Y, Ex, Ey);

axis([xlim(1) xlim(2) ylim(1) ylim(2)]);
axis equal;

xlabel('x (m)');
ylabel('y (m)');
title('Electric Field of Point Charges');

end
