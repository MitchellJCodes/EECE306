function [t,w] = nodes(a,b,N,rule)
% em.quad.nodes
% Return N quadrature nodes and weights on [a,b]
%
% Rules:
% midpoint
% trapz
% simpson
% gauss

    if nargin < 4
        rule = 'midpoint';
    end

    if ~isscalar(N) || N < 1 || N ~= floor(N)
        error('N must be a positive integer');
    end

    rule = lower(rule);

    switch rule

        % Midpoint rule
        case 'midpoint'

            h = (b-a)/N;

            t = a + ((1:N)' - 0.5)*h;
            w = h*ones(N,1);

        % Trapezoidal rule
        % N points including both endpoints
        case 'trapz'

            if N == 1
                t = (a+b)/2;
                w = b-a;
            else
                t = linspace(a,b,N)';
                h = (b-a)/(N-1);

                w = h*ones(N,1);
                w(1) = h/2;
                w(end) = h/2;
            end

        % Composite Simpson rule. N is the requested number of nodes.
        % Odd N uses Simpson 1/3 throughout. Even N uses Simpson 1/3
        % on the first N-3 nodes and Simpson 3/8 on the final 4 nodes,
        % preserving N nodes while retaining fourth-order accuracy.
        case 'simpson'

            if N < 2
                t = (a+b)/2;
                w = b-a;
            else
                t = linspace(a,b,N)';
                h = (b-a)/(N-1);
                w = zeros(N,1);

                if mod(N,2) == 1
                    % Standard composite Simpson 1/3 rule.
                    w(1) = 1;
                    w(end) = 1;
                    w(2:2:N-1) = 4;
                    w(3:2:N-2) = 2;
                    w = (h/3)*w;
                elseif N == 2
                    % Only one interval: trapezoidal fallback.
                    w(:) = h/2;
                else
                    % Simpson 1/3 on N-3 nodes (even number of intervals).
                    n13 = N-3;
                    if n13 >= 3
                        w(1) = h/3;
                        w(2:2:n13-1) = 4*h/3;
                        w(3:2:n13-2) = 2*h/3;
                        w(n13) = h/3;
                    end
                    % Simpson 3/8 on the final three intervals.
                    w(N-3) = w(N-3) + 3*h/8;
                    w(N-2) = 9*h/8;
                    w(N-1) = 9*h/8;
                    w(N) = 3*h/8;
                end
            end

        % Gauss-Legendre rule
        case 'gauss'

            if N == 1

                t = (a+b)/2;
                w = b-a;

            else

                % Jacobi matrix for Gauss-Legendre quadrature
                %
                % beta_k = k/sqrt(4*k^2 - 1)
                k = (1:N-1)';

                beta = k ./ sqrt(4*k.^2 - 1);

                J = diag(beta,1) + diag(beta,-1);

                [V,D] = eig(J);

                x = diag(D);

                [x,idx] = sort(x);
                V = V(:,idx);

                % Weights on [-1,1]
                w0 = 2*(V(1,:).^2)';

                % Map [-1,1] -> [a,b]
                t = (a+b)/2 + (b-a)/2*x;
                w = (b-a)/2*w0;
            end


        otherwise
            error('Unknown quadrature rule: %s',rule);

    end

    % Force required output shapes
    t = t(:);
    w = w(:);

end
