function [pos,w] = line(curve,interval,N)
% em.quad.line
% Quadrature of a parametrized curve
%
% curve - function handle r(t) returning N-by-3 positions
% interval - [a b]
% N - number of quadrature points
%
% Returns:
% pos - N-by-3 positions
% w - N-by-1 physical line-length weights

    a = interval(1);
    b = interval(2);

    % Use midpoint quadrature for the parameter
    [t,wt] = em.quad.nodes(a,b,N,'midpoint');

    % Evaluate curve
    pos = curve(t);
    
    % Make sure positions are N-by-3
    pos = reshape(pos,[],3);

    % Numerical derivative dr/dt
    % Five-point centered difference
    %
    % f'(t) ~= [f(t-2h)-8f(t-h)+8f(t+h)-f(t+2h)]/(12h)
    %
    % h chosen to balance truncation and roundoff
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

    % Speed |dr/dt|
    speed = sqrt(sum(drdt.^2,2));

    % Physical line element
    % ds = |dr/dt| dt
    %
    % abs(wt) makes geometric length positive even if the
    % parameter interval is supplied backwards
    w = abs(wt).*speed;

    w = w(:);

end

