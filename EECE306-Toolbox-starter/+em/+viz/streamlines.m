function h = streamlines(F, xlim, ylim, seeds, opts)
% Visualize streamlines of a vector field in the xy-plane.

% Set default options
if nargin < 5 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'N')
    opts.N = 100;
end

if ~isfield(opts,'z')
    opts.z = 0;
end

% Create observation grid
x = linspace(xlim(1),xlim(2),opts.N);
y = linspace(ylim(1),ylim(2),opts.N);

[X,Y] = meshgrid(x,y);

% Convert grid into observation points
r = [X(:), Y(:), opts.z*ones(numel(X),1)];

% Evaluate vector field in batches
Fv = zeros(size(r));

batchSize = 500;

for k = 1:batchSize:size(r,1)

    idx = k:min(k+batchSize-1,size(r,1));

    Fv(idx,:) = F(r(idx,:));

end

% Extract vector field components
Fx = Fv(:,1);
Fy = Fv(:,2);

% Reshape components to match grid
U = reshape(Fx,size(X));
V = reshape(Fy,size(Y));

% Extract seed coordinates
startX = seeds(:,1);
startY = seeds(:,2);

% Plot streamlines
h = streamline(X,Y,U,V,startX,startY);

% Format plot
axis equal;
axis([xlim ylim]);

xlabel('x (m)');
ylabel('y (m)');
title('Vector Field Streamlines');

grid on;

end