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
        % N is the number of subintervals.
        % N must be even.
        %
        % The rule uses N+1 equally spaced nodes.
        case 'simpson'

            if N == 1
                
                % With one subinterval, Simpson is not possible
                % Fall back to trapezoidal rule
                t = linspace(a,b,2)';
                h = b-a;

                w = [h/2; h/2];

            elseif mod(N,2) ~= 0

                error('Simpson rule requires N to be even');

            else

                h = (b-a)/N;

                t = linspace(a,b,N+1)';

                % Start with the endpoint weights
                w = ones(N+1,1);

                % Interior even-indexed nodes get weight 2
                w(3:2:N-1) = 2;

                % Interior odd-indexed nodes get weight 4
                w(2:2:N) = 4;

                w = (h/3)*w;
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
