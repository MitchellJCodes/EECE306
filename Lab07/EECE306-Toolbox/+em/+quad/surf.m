
function [pos,w] = surf(surface,uinterval,vinterval,Nu,Nv)
% em.quad.surf
% Quadrature over a parametrized surface
%
% surface - r(u,v), returning points in R^3
% uinterval - [u0 u1]
% vinterval - [v0 v1]
% Nu, Nv - number of quadrature points
%
% Returns:
% pos - Nu*Nv by 3
% w - Nu*Nv by 1

    ua = uinterval(1);
    ub = uinterval(2);

    va = vinterval(1);
    vb = vinterval(2);

    % Midpoint quadrature in each parameter
    [u,wu] = em.quad.nodes(ua,ub,Nu,'midpoint');
    [v,wv] = em.quad.nodes(va,vb,Nv,'midpoint');
    
    % Tensor-product grid
    [U,V] = ndgrid(u,v);

    U = U(:);
    V = V(:);

    % Corresponding tensor-product parameter weights
    [WU,WV] = ndgrid(wu,wv);

    wt = abs(WU(:).*WV(:));

    % Surface positions
    pos = surface(U,V);
    pos = reshape(pos,[],3);

    % Numerical partial derivatives.
    hu = eps^(1/5) .* (1 + abs(U));
    hv = eps^(1/5) .* (1 + abs(V));

    r_um2 = surface(U - 2*hu,V);
    r_um1 = surface(U - hu,V);
    r_up1 = surface(U + hu,V);
    r_up2 = surface(U + 2*hu,V);

    r_vm2 = surface(U,V - 2*hv);
    r_vm1 = surface(U,V - hv);
    r_vp1 = surface(U,V + hv);
    r_vp2 = surface(U,V + 2*hv);

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
