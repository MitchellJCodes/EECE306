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

        % Simpson rule
        %
        % For odd N:
        % standard composite Simpson 1/3.
        %
        % For even N:
        % Simpson 1/3 on the first N-3 point
        % Simpson 3/8 on the final 4 points
        % This preserves fourth-order convergence
        case 'simpson'

            if N == 1
                t = (a+b)/2;
                w = b-a;

            elseif N == 2

                % With only two points, Simpson is not possible
                % Fall back to trapezoidal rule
                t = linspace(a,b,N)';
                h = b-a;

                w = [h/2; h/2];

            else

                t = linspace(a,b,N)';
                h = (b-a)/(N-1);

                w = zeros(N,1);

                if mod(N,2) == 1

                    % Odd number of points -> even number of intervals
                    % Standard Simpson 1/3
                    w(1) = 1;
                    w(N) = 1;

                    w(2:2:N-1) = 4;
                    w(3:2:N-2) = 2;

                    w = (h/3)*w;

                else

                    % Even number of points -> odd number of intervals
                    %
                    % Simpson 1/3 over the first N-3 points
                    % Simpson 3/8 over the last 4 points

                    m = N-3;

                    % First section: points 1:m
                    w(1) = w(1) + 1;
                    w(m) = w(m) + 1;

                    if m >= 3
                        w(2:2:m-1) = w(2:2:m-1) + 4;
                        w(3:2:m-2) = w(3:2:m-2) + 2;
                    end

                    w = (h/3)*w;

                    % Last three intervals: points m:N
                    %
                    % 3/8 weights = [1 3 3 1] * 3h/8
                    w(m:N) = w(m:N) + ...
                             (3*h/8)*[1; 3; 3; 1];
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
