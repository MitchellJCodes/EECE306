function [p, h] = convergence(errFcn, Nlist, opts)
%CONVERGENCE Convergence Test
%   [P, H] = EM.TEST.CONVERGENCE(ERRFCN, NLIST, OPTS) returns the observed
%   order p and a log-log plot.
%
%   Example: errFcn = @(N) 1./N.^2;
%   Nlist = [10 20 40 80 160];
%   opts.marker = 'o-';
%   [p, h] = em.test.convergence(errFcn, Nlist, opts); % returns p = 2.0000
%   and the linear error plot

if nargin < 3
    opts = struct();
end

% Error for each N
E = zeros(size(Nlist));
for k = 1:numel(Nlist)
    E(k) = errFcn(Nlist(k));
end

% Log-log regression with flat vectors
x = log(Nlist(:));
y = log(E(:));

% Slope and intercepts from best fit polynomial
c = polyfit(x, y, 1);
slope = c(1);

% Observed order
p = -slope;

% Plot marker default
if ~isfield(opts, 'marker')
    opts.marker = 'o-';
end

% Plot creation
h = figure;
loglog(Nlist, E, opts.marker, 'LineWidth', 1.5);
grid on
xlabel('N')
ylabel('error')

if isfield(opts, 'title')
    title(opts.title)
else
    title(sprintf('Observed order p = %3f', p));
end

end

