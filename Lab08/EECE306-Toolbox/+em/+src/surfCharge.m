function s = surfCharge(rhoS,surf,uspan,vspan,Nu,Nv,rule)
% em.src.surfCharge
% Construct a discretized surface-charge source
%
% rhoS - scalar surface charge density or function handle rhoS(r)
% surf - parametrized surface r(u,v)
% uspan - [u0 u1]
% vspan - [v0 v1]
% Nu,Nv - number of quadrature points
% rule - quadrature rule
%
% Source fields:
% s.type
% s.pos
% s.w
% s.q
% s.Idl
% s.label

    if nargin < 7 || isempty(rule)
        rule = 'midpoint';
    end

    % Get positions and physical surface-area weights
    [pos,w] = em.quad.surf(surf,uspan,vspan,Nu,Nv,rule);

    % Evaluate surface charge density
    if isa(rhoS,'function_handle')
        q = rhoS(pos);
    else
        q = rhoS*ones(size(w));
    end

    q = q(:);
    w = w(:);

    if numel(q) ~= Nu*Nv
        error('rhoS must return one value per quadrature point');
    end

    % Build source structure
    s.type = 'charge';
    s.pos = pos;
    s.w = w;
    s.q = q;
    s.Idl = [];
    s.label = '';

end
