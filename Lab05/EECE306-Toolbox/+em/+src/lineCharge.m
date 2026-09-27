function s = lineCharge(rhoL,curve,tspan,N,rule)
% em.src.lineCharge
% Construct a discretized line-charge source
%
% rhoL - scalar charge density or function handle rhoL(r)
% curve - parametrized curve r(t)
% tspan - [t0 t1]
% N - number of quadrature points
% rule - quadrature rule
%
% Source fields:
% s.type
% s.pos
% s.w
% s.q
% s.Idl
% s.label

    if nargin < 5 || isempty(rule)
        rule = 'midpoint';
    end

    % Get positions and physical line-element weights
    [pos,w] = em.quad.line(curve,tspan,N,rule);

    % Evaluate charge density
    if isa(rhoL,'function_handle')
        q = rhoL(pos);
    else
        q = rhoL*ones(size(w));
    end

    q = q(:);
    w = w(:);

    if numel(q) ~= N
        error('rhoL must return one value per quadrature point');
    end

    % Build source structure
    s.type = 'charge';
    s.pos = pos;
    s.w = w;
    s.q = q;
    s.Idl = [];
    s.label = '';

end
