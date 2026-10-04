function [t,w] = nodes(a,b,N,rule)
% Defining midpoint as the default rule
if nargin < 4
    rule = 'midpoint';
end

% t explains where to sample
% w explains how much each sample counts

switch rule

    case 'midpoint'
        h = (b-a)/N;
        
        t = (a + h/2 : h : b - h/2).'; % first midpoint, moves by an increments of h, and stops at the final midpoint

        w = h * ones(N,1);

    case 'trapz'
        if N == 1
            t = (a+b)/2;
            w = b-a;
        else
            t = linspace(a,b,N).';
            h = (b-a)/(N-1);
            w = h * ones(N,1);
        
            w(1) = h/2;
            w(end) = h/2;
        end
    case 'simpson'
         if mod(N,2) ~= 0
            error('N must be even for Simpson rule');
        end

        % N intervals means N+1 nodes
        h = (b-a)/N;
        t = linspace(a,b,N+1).';

        % Simpson weights: 1 4 2 4 2 ... 4 1
        w = ones(N+1,1);
        w(2:2:end-1) = 4;
        w(3:2:end-2) = 2;

        w = (h/3)*w;

    case 'gauss'
        beta = (1:N-1) ./ sqrt(4*(1:N-1).^2 -1);

        J = diag(beta,1) + diag(beta,-1);

        [V,D] = eig(J);

        x = diag(D);

        [x, idx] = sort(x);

        V = V(:, idx);
        
        w_ref = 2 * (V(1,:).^2).';

        % map from  [-1,1] to [a,b]
        t = (a+b)/2 + (b-a)/2 * x;
        w = (b-a)/2 * w_ref;
    otherwise
        error('Unknown rule specified');
end

end

