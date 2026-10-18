function h = quiver2(F, xlim, ylim, N, opts)
%QUIVER2 Electric Field Plot Function
%   H = EM.VIZ.QUIVER2(F, XLIM, YLIM, N, OPTS) returns the vectorized plot 
%   of an electric field given a function handle to sample.
%
%   Example:
%   figure;
%   em.viz.quiver2(@(r) em.field.E(s2, r), [-1 1], [-1 1], 15, ...
%   struct('normalize', false, 'title', 'Two equal charges, raw plot'));
%   % returns an arrow plot of a vector field on a plane, with the title 
%     'Two equal charges, raw plot'
%
%   NOTE: if false is instead set to true, it will normalize the plot and
%   makes the field directions easier to read. However, the relative
%   magnitudes of the field lines will be lost.

if nargin < 5 || isempty(opts)
    opts = struct();
end
if ~isfield(opts,'normalize') || isempty(opts.normalize)
    opts.normalize = false;
end
if ~isfield(opts,'title') || isempty(opts.title)
    opts.title = 'Electric Field';
end
if ~isfield(opts,'newFigure') || isempty(opts.newFigure)
    opts.newFigure = true;
end
if ~isa(F,'function_handle')
    error('F must be a function handle');
end
if ~isscalar(N) || N < 2 || N ~= floor(N)
    error('N must be an integer greater than or equal to 2');
end
if numel(xlim) ~= 2 || numel(ylim) ~= 2
    error('xlim and ylim must each contain two values');
end

% Create N evenly spaced x and y coordinates within the plot limits.
x = linspace(xlim(1), xlim(2), N);
y = linspace(ylim(1), ylim(2), N);
[X,Y] = meshgrid(x, y);
r = [X(:), Y(:), zeros(numel(X),1)];
E = F(r);
if ~isequal(size(E),size(r))
    error('F must return an N-by-3 array matching its N-by-3 input');
end
Ex = reshape(E(:,1), size(X));
Ey = reshape(E(:,2), size(Y));

if opts.normalize
    M = hypot(Ex,Ey);
    M(M == 0) = 1; % avoid division by zero while preserving zero arrows
    Ex = Ex ./ M;
    Ey = Ey ./ M;
end

if opts.newFigure
    figure;
else
    % Preserve existing contour/field graphics when explicitly overlaying.
    hold on;
end
quiver(X, Y, Ex, Ey);
axis([xlim(1) xlim(2) ylim(1) ylim(2)]);
axis equal;
xlabel('x (m)');
ylabel('y (m)');
title(opts.title);
if ~opts.newFigure
    hold off;
end
end