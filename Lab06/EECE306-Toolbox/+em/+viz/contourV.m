function contourV(V, xlim, ylim, N, opts)
% CONTOURV Plot contours of a scalar electric potential.
%
%   em.viz.contourV(V,xlim,ylim,N,opts)
%
%   V     - scalar field function handle
%   xlim  - [xmin xmax] plot limits
%   ylim  - [ymin ymax] plot limits
%   N     - number of points in each coordinate direction
%   opts  - options structure
%
%   The scalar field is evaluated in the z = 0 plane and displayed
%   using a filled contour plot.
%
%   opts.title - title for the plot

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
Vval = V(r);

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
