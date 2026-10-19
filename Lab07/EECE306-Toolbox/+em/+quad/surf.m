
function [pos,w] = surf(surface,uinterval,vinterval,Nu,Nv,rule)
% em.quad.surf
% Quadrature over a parametrized surface
%
% surface - r(u,v), returning points in R^3
% uinterval - [u0 u1]
% vinterval - [v0 v1]
% Nu, Nv - number of quadrature points
% rule - quadrature rule (optional; defaults to midpoint)
%
% Returns:
% pos - Nu*Nv by 3
% w - Nu*Nv by 1

    if nargin < 6 || isempty(rule)
        rule = 'midpoint';
    end

    ua = uinterval(1);
    ub = uinterval(2);

    va = vinterval(1);
    vb = vinterval(2);

    % Midpoint quadrature in each parameter
    [u,wu] = em.quad.nodes(ua,ub,Nu,rule);
    [v,wv] = em.quad.nodes(va,vb,Nv,rule);
    
    % Tensor-product grid
    [U,V] = ndgrid(u,v);

    U = U(:);
    V = V(:);

    % Corresponding tensor-product parameter weights
    [WU,WV] = ndgrid(wu,wv);

    wt = abs(WU(:).*WV(:));

    % Surface positions
    pos = evalSurface(surface,U,V);

    % Numerical partial derivatives.
    hu = eps^(1/5) .* (1 + abs(U));
    hv = eps^(1/5) .* (1 + abs(V));

    r_um2 = evalSurface(surface,U - 2*hu,V);
    r_um1 = evalSurface(surface,U - hu,V);
    r_up1 = evalSurface(surface,U + hu,V);
    r_up2 = evalSurface(surface,U + 2*hu,V);

    r_vm2 = evalSurface(surface,U,V - 2*hv);
    r_vm1 = evalSurface(surface,U,V - hv);
    r_vp1 = evalSurface(surface,U,V + hv);
    r_vp2 = evalSurface(surface,U,V + 2*hv);

    r_um2 = reshape(r_um2,[],3);
    r_um1 = reshape(r_um1,[],3);
    r_up1 = reshape(r_up1,[],3);
    r_up2 = reshape(r_up2,[],3);

    r_vm2 = reshape(r_vm2,[],3);
    r_vm1 = reshape(r_vm1,[],3);
    r_vp1 = reshape(r_vp1,[],3);
    r_vp2 = reshape(r_vp2,[],3);

    ru = (r_um2 - 8*r_um1 + 8*r_up1 - r_up2) ./ (12*hu);
    rv = (r_vm2 - 8*r_vm1 + 8*r_vp1 - r_vp2) ./ (12*hv);

    % Surface Jacobian
    % dS = |ru x rv| du dv
    crossprod = zeros(size(ru));

    crossprod(:,1) = ru(:,2).*rv(:,3) - ru(:,3).*rv(:,2);
    crossprod(:,2) = ru(:,3).*rv(:,1) - ru(:,1).*rv(:,3);
    crossprod(:,3) = ru(:,1).*rv(:,2) - ru(:,2).*rv(:,1);

    jac = sqrt(sum(crossprod.^2,2));

    w = wt .* jac;

    w = w(:);

end

function r = evalSurface(fun,u,v)
% Support both vectorized element-wise handles and scalar-only handles.
try
    r = fun(u,v);
    if numel(r) == 3*numel(u)
        r = reshape(r,[],3);
        return
    end
catch
    % Retry pointwise below for parameterizations written with * or ^.
end
r = zeros(numel(u),3);
for ii = 1:numel(u)
    ri = fun(u(ii),v(ii));
    if numel(ri) ~= 3
        error('surface must return one 1-by-3 position per parameter pair');
    end
    r(ii,:) = reshape(ri,1,3);
end
end
