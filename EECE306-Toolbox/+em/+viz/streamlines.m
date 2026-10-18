function h = streamlines(F, xlim, ylim, seeds, opts)
%STREAMLINES Field Line Sketches
%   H = STREAMLINES(F, XLIM, YLIM, SEEDS, OPTS) sketches the field lines
%   from the seed points.
%   
%   Example: FE = @(r) em.field.E(s, r);
%   em.viz.streamlines(FE, [-2 2], [-2 2], [cosd(0:45:315)' sind(0:45:315)' zeros(8,1)].*0.2, struct());
%   title('Equipotentials with field lines overlaid') % returns the quiver2
%   plot with field lines over it.

if isempty(fieldnames(opts))
    opts.normalize = false;
    opts.StepSize =  0.5;
    opts.MaxSteps = 500;
end

hold on

h = gobjects(size(seeds,1),1);

for n = 1:size(seeds,1)

    seed = seeds(n,:);

    % Forward trace
    rf = zeros(opts.MaxSteps,3);
    rf(1,:) = seed;
    nf = 1;

    for k = 1:opts.MaxSteps-1

        v = F(rf(k,:));

        if norm(v) < eps
            break
        end

        v = v / norm(v);

        rf(k+1,:) = rf(k,:) + opts.StepSize*v;

        x = rf(k+1,1);
        y = rf(k+1,2);

        nf = k + 1;

        if x < xlim(1) || x > xlim(2) || ...
           y < ylim(1) || y > ylim(2)
            break
        end

    end

    rf = rf(1:nf,:);

    % Backward trace
    rb = zeros(opts.MaxSteps,3);
    rb(1,:) = seed;
    nb = 1;

    for k = 1:opts.MaxSteps-1

        v = -F(rb(k,:));

        if norm(v) < eps
            break
        end

        v = v / norm(v);

        rb(k+1,:) = rb(k,:) + opts.StepSize*v;

        x = rb(k+1,1);
        y = rb(k+1,2);

        nb = k + 1;

        if x < xlim(1) || x > xlim(2) || ...
           y < ylim(1) || y > ylim(2)
            break
        end

    end

    rb = rb(1:nb,:);

    % Join forward and backward traces
    rb = flipud(rb);

    r = [rb(1:end-1,:); rf];

    h(n) = plot(r(:,1), r(:,2), 'k-');

end

end