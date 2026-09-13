function [t,w] = nodes(a,b,N,rule)
%NODES 1-D quadrature nodes and weights on [a,b].

    if nargin < 4
        rule = 'midpoint';
    end

    switch lower(rule)

        case 'midpoint'

            h = (b-a)/N;

            t = a + ((1:N)' - 0.5)*h;
            w = h*ones(N,1);

        otherwise

            error('Unknown quadrature rule: %s', rule);

    end

end
