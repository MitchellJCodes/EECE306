function h = contourV(F, xlim, ylim, N, opts)
%CONTOURV Contour Plot
%   H = EM.VIZ.COUNTOURV(F, XLIM, YLIM, N, OPTS) returns a contour plot of
%   a scalar field.
%
%   Example: s  = em.src.pointCharge(Q, [0 0 0]);
%   figure;
%   em.viz.contourV(@(r) em.field.V(s, r), [-3 3], [-3 3], 60, struct('title', 'Equipotentials of a point charge'));

if nargin < 5 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'title')
    opts.title = 'Electric Potential';
end

% Create N evenly spaced x and y coordinates
x = linspace(xlim(1),xlim(2),N);
y = linspace(ylim(1),ylim(2),N);

% Create a 2-D grid
[X,Y] = meshgrid(x,y);

% Convert grid to N-by-3 observation points in the z = 0 plane
r = [X(:), Y(:), zeros(numel(X),1)];

% Evaluate scalar potential
Vval = F(r);

% Reshape potential values to match the grid
Vval = reshape(Vval,size(X));

% Plot filled equipotential contours
figure;
contourf(X,Y,Vval,20);
colorbar;

axis([xlim(1) xlim(2) ylim(1) ylim(2)]);
axis equal;

xlabel('x (m)');
ylabel('y (m)');
title(opts.title);

end

