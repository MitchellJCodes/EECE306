function h = contourV(f, xlim, ylim, N, opts)

if nargin < 5
    opts = struct();
end

% Create N evenly spaced x and y coordinates
x = linspace(xlim(1),xlim(2),N);
y = linspace(ylim(1),ylim(2),N);

% Create an N x N grid
[X,Y] = meshgrid(x,y);

% Convert grid into N^2 x 3 observation points
r = [X(:),Y(:),zeros(numel(X),1)];

% Calculate potential at each observation point
V = f(r);

% Reshape potential values to match the grid
Vgrid = reshape(V,size(X));

% Create filled contour plot
figure;
h = contourf(X,Y,Vgrid,20);

colorbar;
hold on;

axis([xlim(1) xlim(2) ylim(1) ylim(2)]);
axis equal;

xlabel('x (m)');
ylabel('y (m)');
title('Electric Potential Contours');

end