function [t, w] = nodes(a, b, N, rule)
%NODES Nodes and Weights Engine
%   [T, W] = EM.QUAD.NODES(A, B, N, RULE) returns a Nx1 column vector of 
%   nodes (T) and an Nx1 column vector of weights (W).
%   
%   Example: [t, w] = em.quad.nodes(0, 1, 8, rule) % returns two vectors:
%   t = [0.1250,        w = [0.2500,
%        0.3750,             0.2500,
%        0.6250,             0.2500,
%        0.8750]             0.2500]
%
%   NOTE 1: The Midpoint rule is the default rule if no rule is given.
%   NOTE 2: The Trapezoid rule and Simpson rule take N+1 sample points.
%   NOTE 3: The standard Simpson rule requires an even, positive number of 
%   intervals because it uses pairs of subintervals.


if nargin < 4 % Default rule = midpoint rule when rule is omitted
    rule = 'midpoint';
end

if N <= 0 % Guarding the singularity
    error('ERROR: Invalid number of intervals')
end

rule = lower(rule);

switch rule % Case statement for each rule

    case 'midpoint' % Midpoint Rule
        h = (b - a)/N; % Width of interval
        t = a + ((1:N)' - 0.5) * h; % Nodes from first to last
        w = h * ones(N, 1); % Vector of weights (single weight = h)

    case 'trapz' % Trapezoid Rule
        if N == 1 % For one subinterval...
            t = (a+b)/2; % ...the sample point is the midpoint...
            w = b-a; % ...and weight is the length of the range.
        else % Case with many subintervals
            h = (b-a)/(N-1); % Interval width; (N-1) takes care of N+1 node
            t = linspace(a, b, N).'; % Evenly spaced nodes
            w = h * ones (N, 1); % Interior weights
            w(1) = h/2; % First weight is width / 2
            w(end) = h/2; % Last weight is also width / 2
        end

    case 'simpson' %Composite Simpson rule
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
      

    case 'gauss' % Gauss Rule - Golub-Welsch approach
        if N == 1
            t = (a+b)/2;
            w = b-a;

        else
        k = (1 : N-1)'; % Column vector of integers from 1 to N-1
        beta = k ./ sqrt(4*k.^2 - 1); % Vector of coefficients
        J = diag(beta, 1) + diag(beta, -1); % Jacobi matrix
        % D = diagonal matrix of eigenvalues
        % V = matrix with columns eigenvector 
        [V, D] = eig(J);
        x = diag(D); % Nodes on [-1, 1] (will be irrational)
        [x, idx] = sort(x); % Sorts nodes
        V = V(:, idx); % Sorts eigenvectors
        % V(1,:) = first components of eigenvectors
        w_ref = 2*(V(1,:).^2)'; % Weights on [-1, 1]
        t = ((a+b) + x*(b-a))/2; % Mapping nodes to [a, b]
        w = w_ref*(b-a)/2; % Mapping weights to [a, b]
        end

    otherwise
        error('ERROR: Invalid input. Rule not found.'\n)
        printf('Valid rules: midpoint, trapz, simpson, gauss.')

end

% Force required output shapes
t = t(:);
w = w(:);

end

