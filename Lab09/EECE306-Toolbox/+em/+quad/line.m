function [pos,w,drdt] = line(curve,interval,N,rule)
% em.quad.line
% Quadrature of a parametrized curve
%
% curve - function handle r(t) returning N-by-3 positions
% interval - [a b]
% N - number of quadrature points
% rule - quadrature rule (optional; defaults to midpoint)
%
% Returns:
% pos - N-by-3 positions
% w - N-by-1 physical line-length weights
% drdt - N-by-3 parameter tangent dr/dt (optional third output)

    if nargin < 4 || isempty(rule)
        rule = 'midpoint';
    end

    a = interval(1);
    b = interval(2);

    % Midpoint
    [t,wt] = em.quad.nodes(a,b,N,rule);

    % Curve
    pos = curve(t);
    
    % Positions N-by-3
    pos = reshape(pos,[],3);

    % Numerical derivative dr/dt
    h = eps^(1/5) .* (1 + abs(t));

    rm2 = curve(t - 2*h);
    rm1 = curve(t - h);
    rp1 = curve(t + h);
    rp2 = curve(t + 2*h);

    rm2 = reshape(rm2,[],3);
    rm1 = reshape(rm1,[],3);
    rp1 = reshape(rp1,[],3);
    rp2 = reshape(rp2,[],3);

    drdt = (rm2 - 8*rm1 + 8*rp1 - rp2) ./ (12*h);

    speed = sqrt(sum(drdt.^2,2));

    w = abs(wt).*speed;

    w = w(:);

end

