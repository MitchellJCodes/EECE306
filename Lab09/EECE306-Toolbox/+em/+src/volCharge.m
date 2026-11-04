function s = volCharge(rhoV,vol,uspan,vspan,wspan,Nu,Nv,Nw,rule)
% em.src.volCharge
% Construct a discretized volume-charge source
%
% rhoV - scalar volume charge density or function handle rhoV(r)
% vol - parametrized volume r(u,v,w)
% uspan - [u0 u1]
% vspan - [v0 v1]
% wspan - [w0 w1]
% Nu,Nv,Nw - number of quadrature points
% rule - quadrature rule
%
% Source fields:
% s.type
% s.pos
% s.w
% s.q
% s.Idl
% s.label

    if nargin < 9 || isempty(rule)
        rule = 'midpoint';
    end

    % Get positions and physical volume weights
    [pos,w] = em.quad.vol( ...
        vol,uspan,vspan,wspan,Nu,Nv,Nw,rule);

    % Evaluate volume charge density
    if isa(rhoV,'function_handle')
        q = rhoV(pos);
    else
        q = rhoV*ones(size(w));
    end

    q = q(:);
    w = w(:);

    if numel(q) ~= Nu*Nv*Nw
        error('rhoV must return one value per quadrature point');
    end

    % Build source structure
    s.type = 'charge';
    s.pos = pos;
    s.w = w;
    s.q = q;
    s.Idl = [];
    s.label = '';

end

