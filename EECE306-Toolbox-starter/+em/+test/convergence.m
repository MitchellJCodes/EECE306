function p =  convergence(errFcn, Nlist, opts)
% errFcn  = a function that calculates the error for a given N
% Nlist   = different N values you want to test
% opts    = options, such as whether to make a plot

errors = zeros(size(Nlist));

for k = 1:length(Nlist)
    errors(k) = errFcn(Nlist(k));
end

fit = polyfit(log(Nlist),log(errors),1);

p = -fit(1);

if opts.plot
    figure;
    loglog(Nlist,errors,'o-');
    xlabel('N');
    ylabel('Error');
    title(sprintf('Observed order p = %.2f',p));
    grid on;
end

end