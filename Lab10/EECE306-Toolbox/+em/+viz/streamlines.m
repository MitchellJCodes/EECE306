function lines = streamlines(F, xlim, ylim, seeds, opts)
%STREAMLINES Trace and plot planar vector-field streamlines.
%
%   em.viz.streamlines(F,xlim,ylim,seeds,opts)
%   lines = em.viz.streamlines(...)
%
%   F      - vector field function handle accepting N-by-3 points and
%            returning N-by-3 vectors. The x-y components are traced in z=0.
%   xlim   - [xmin xmax] plot limits
%   ylim   - [ymin ymax] plot limits
%   seeds  - N-by-2 array of starting points in the x-y plane
%   opts   - optional struct:
%              stepSize       fixed arc-length step (default min(span)/250)
%              maxSteps       maximum steps per line (default 5000)
%              maxArcLength   maximum traced length (default 4*domain diagonal)
%              closeTolerance distance for detecting return to seed
%                             (default 2.5*stepSize)
%              minSteps       minimum steps before closure test (default 20)
%              minField       field magnitude below which tracing stops
%                             (default 1e-12)
%              lineWidth      plotted line width (default 1.2)
%              title          axes title (default 'Magnetic Field Lines')
%              newFigure      create a new figure (default false)
%
%   lines  - cell array containing the traced N-by-2 point arrays.
%
%   A normalized field direction is integrated with fixed-step RK4. The
%   field is normalized so that stepSize controls geometric spacing rather
%   than the magnitude of F. Tracing stops at the plot boundary, a field
%   null, the maximum length, or when a line returns to its seed.

if nargin < 5 || isempty(opts)
    opts = struct();
end
if ~isa(F, 'function_handle')
    error('em.viz.streamlines:InvalidField', ...
        'F must be a function handle.');
end
if ~isstruct(opts) || numel(opts) ~= 1
    error('em.viz.streamlines:InvalidOptions', ...
        'opts must be a scalar struct.');
end
if ~isnumeric(xlim) || numel(xlim) ~= 2 || any(~isfinite(xlim(:))) || ...
        xlim(2) <= xlim(1)
    error('em.viz.streamlines:InvalidXLim', ...
        'xlim must be a finite increasing two-element vector.');
end
if ~isnumeric(ylim) || numel(ylim) ~= 2 || any(~isfinite(ylim(:))) || ...
        ylim(2) <= ylim(1)
    error('em.viz.streamlines:InvalidYLim', ...
        'ylim must be a finite increasing two-element vector.');
end
if ~isnumeric(seeds) || size(seeds,2) ~= 2 || isempty(seeds) || ...
        any(~isfinite(seeds(:)))
    error('em.viz.streamlines:InvalidSeeds', ...
        'seeds must be a nonempty finite N-by-2 array.');
end
if any(seeds(:,1) < xlim(1) | seeds(:,1) > xlim(2) | ...
       seeds(:,2) < ylim(1) | seeds(:,2) > ylim(2))
    error('em.viz.streamlines:SeedOutsideDomain', ...
        'Every seed must lie inside the plotting limits.');
end

span = min([diff(xlim), diff(ylim)]);
opts = setDefault(opts, 'stepSize', span/250);
opts = setDefault(opts, 'maxSteps', 5000);
opts = setDefault(opts, 'maxArcLength', ...
    4*hypot(diff(xlim), diff(ylim)));
opts = setDefault(opts, 'closeTolerance', 2.5*opts.stepSize);
opts = setDefault(opts, 'minSteps', 20);
opts = setDefault(opts, 'minField', 1e-12);
opts = setDefault(opts, 'lineWidth', 1.2);
opts = setDefault(opts, 'title', 'Magnetic Field Lines');
opts = setDefault(opts, 'newFigure', false);

if ~isscalar(opts.stepSize) || ~isfinite(opts.stepSize) || opts.stepSize <= 0 || ...
   ~isscalar(opts.maxSteps) || opts.maxSteps < 1 || ...
        opts.maxSteps ~= floor(opts.maxSteps) || ...
   ~isscalar(opts.maxArcLength) || ~isfinite(opts.maxArcLength) || ...
        opts.maxArcLength <= 0 || ...
   ~isscalar(opts.closeTolerance) || ~isfinite(opts.closeTolerance) || ...
        opts.closeTolerance <= 0 || ...
   ~isscalar(opts.minSteps) || opts.minSteps < 1 || ...
        opts.minSteps ~= floor(opts.minSteps) || ...
   ~isscalar(opts.minField) || ~isfinite(opts.minField) || opts.minField < 0 || ...
   ~isscalar(opts.lineWidth) || ~isfinite(opts.lineWidth) || opts.lineWidth <= 0
    error('em.viz.streamlines:InvalidOptions', ...
        'Numeric options must be finite and have valid positive values.');
end
if ~ischar(opts.title)
    error('em.viz.streamlines:InvalidTitle', 'opts.title must be a character vector.');
end
if ~isscalar(opts.newFigure)
    error('em.viz.streamlines:InvalidFigureOption', ...
        'opts.newFigure must be a logical scalar.');
end

% Trace each seed independently.
lines = cell(size(seeds,1),1);
for k = 1:size(seeds,1)
    lines{k} = traceOne(F, seeds(k,:), xlim, ylim, opts);
end

% Plot into the current axes unless a new figure was explicitly requested.
if opts.newFigure
    figure;
end
wasHolding = ishold;
hold on;
for k = 1:numel(lines)
    if size(lines{k},1) >= 2
        plot(lines{k}(:,1), lines{k}(:,2), 'LineWidth', opts.lineWidth);
    end
end
axis([xlim(1) xlim(2) ylim(1) ylim(2)]);
axis equal;
xlabel('x (m)');
ylabel('y (m)');
title(opts.title);
if ~wasHolding
    hold off;
end
end

function line = traceOne(F, seed, xlim, ylim, opts)
step = opts.stepSize;
points = zeros(opts.maxSteps+2, 2);
points(1,:) = seed;
count = 1;
arcLength = 0;

for k = 1:opts.maxSteps
    p = points(count,:);
    k1 = directionAt(F, p, opts.minField);
    if isempty(k1), break; end
    k2 = directionAt(F, p + 0.5*step*k1, opts.minField);
    if isempty(k2), break; end
    k3 = directionAt(F, p + 0.5*step*k2, opts.minField);
    if isempty(k3), break; end
    k4 = directionAt(F, p + step*k3, opts.minField);
    if isempty(k4), break; end

    next = p + (step/6)*(k1 + 2*k2 + 2*k3 + k4);
    if next(1) < xlim(1) || next(1) > xlim(2) || ...
       next(2) < ylim(1) || next(2) > ylim(2)
        break;
    end

    count = count + 1;
    points(count,:) = next;
    arcLength = arcLength + step;

    % A closed field line returns to its seed after traversing a loop.
    if k >= opts.minSteps && norm(next-seed) <= opts.closeTolerance
        count = count + 1;
        points(count,:) = seed;
        break;
    end
    if arcLength >= opts.maxArcLength
        break;
    end
end

line = points(1:count,:);
end

function direction = directionAt(F, p, minField)
value = F([p(1) p(2) 0]);
if ~isnumeric(value) || ~isreal(value) || numel(value) ~= 3 || ...
        any(~isfinite(value(:)))
    error('em.viz.streamlines:InvalidFieldOutput', ...
        'F must return one finite real three-component vector.');
end
value = reshape(value,1,3);
v = value(1:2);
m = hypot(v(1),v(2));
if m <= minField
    direction = [];
else
    direction = v./m;
end
end

function opts = setDefault(opts, name, value)
if ~isfield(opts,name) || isempty(opts.(name))
    opts.(name) = value;
end
end
