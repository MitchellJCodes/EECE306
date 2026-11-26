function [p,h] = convergence(errFcn,Nlist,opts)
%CONVERGENCE Estimate observed convergence order from an error function.
%   [P,H] = EM.TEST.CONVERGENCE(ERRFCN,NLIST,OPTS) evaluates ERRFCN(N)
%   for each N in NLIST, fits a line to log(error) versus log(N), and
%   returns the observed order P (the negative fitted slope) and figure
%   handle H. ERRFCN must return one positive scalar error for each N.
%
%   OPTS is optional and may contain a char field 'title'.

if nargin < 2 || isempty(Nlist)
    error('Nlist must contain at least two grid sizes');
end
if nargin < 3 || isempty(opts)
    opts = struct();
end
Nlist = Nlist(:);
if numel(Nlist) < 2 || any(~isfinite(Nlist)) || any(Nlist <= 0)
    error('Nlist must contain at least two positive finite values');
end
if ~isa(errFcn,'function_handle')
    error('errFcn must be a function handle');
end

err = zeros(size(Nlist));
for k = 1:numel(Nlist)
    ek = errFcn(Nlist(k));
    if ~isscalar(ek) || ~isfinite(ek) || ek <= 0
        error('errFcn must return a positive finite scalar error');
    end
    err(k) = ek;
end

fit = polyfit(log(Nlist),log(err),1);
p = -fit(1);
h = figure;
loglog(Nlist,err,'o-', 'LineWidth',1.2);
grid on;
xlabel('N');
ylabel('Error');
if isfield(opts,'title') && ~isempty(opts.title)
    title(opts.title);
else
    title(sprintf('Observed convergence order p = %.3f',p));
end
end
