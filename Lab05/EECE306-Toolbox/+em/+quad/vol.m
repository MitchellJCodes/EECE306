function [pos,w] = vol(volume,uinterval,vinterval,qinterval,Nu,Nv,Nq)
% em.quad.vol
% Quadrature over a parametrized volume
%
% volume - r(u,v,q), returning points in R^3
% uinterval - [u0 u1]
% vinterval - [v0 v1]
% qinterval - [q0 q1]
% Nu,Nv,Nq - number of quadrature points
%
% Returns:
% pos - Nu*Nv*Nq by 3
% w - Nu*Nv*Nq by 1

    ua = uinterval(1);
    ub = uinterval(2);

    va = vinterval(1);
    vb = vinterval(2);

    qa = qinterval(1);
    qb = qinterval(2);

    % Midpoint quadrature in each parameter
    [u,wu] = em.quad.nodes(ua,ub,Nu,'midpoint');
    [v,wv] = em.quad.nodes(va,vb,Nv,'midpoint');
    [q,wq] = em.quad.nodes(qa,qb,Nq,'midpoint');

    % Tensor-product grid
    [U,V,Q] = ndgrid(u,v,q);

    U = U(:);
    V = V(:);
    Q = Q(:);

    % Tensor-product parameter weights
    [WU,WV,WQ] = ndgrid(wu,wv,wq);

    wt = abs(WU(:).*WV(:).*WQ(:));

    % Positions
    pos = volume(U,V,Q);
    pos = reshape(pos,[],3);

    % Numerical derivatives.
    hu = eps^(1/5) .* (1 + abs(U));
    hv = eps^(1/5) .* (1 + abs(V));
    hq = eps^(1/5) .* (1 + abs(Q));

    % Derivative with respect to u
    r_um2 = volume(U-2*hu,V,Q);
    r_um1 = volume(U-hu,V,Q);
    r_up1 = volume(U+hu,V,Q);
    r_up2 = volume(U+2*hu,V,Q);

    % Derivative with respect to v
    r_vm2 = volume(U,V-2*hv,Q);
    r_vm1 = volume(U,V-hv,Q);
    r_vp1 = volume(U,V+hv,Q);
    r_vp2 = volume(U,V+2*hv,Q);

    % Derivative with respect to q
    r_qm2 = volume(U,V,Q-2*hq);
    r_qm1 = volume(U,V,Q-hq);
    r_qp1 = volume(U,V,Q+hq);
    r_qp2 = volume(U,V,Q+2*hq);

    r_um2 = reshape(r_um2,[],3);
    r_um1 = reshape(r_um1,[],3);
    r_up1 = reshape(r_up1,[],3);
    r_up2 = reshape(r_up2,[],3);

    r_vm2 = reshape(r_vm2,[],3);
    r_vm1 = reshape(r_vm1,[],3);
    r_vp1 = reshape(r_vp1,[],3);
    r_vp2 = reshape(r_vp2,[],3);

    r_qm2 = reshape(r_qm2,[],3);
    r_qm1 = reshape(r_qm1,[],3);
    r_qp1 = reshape(r_qp1,[],3);
    r_qp2 = reshape(r_qp2,[],3);

    ru = (r_um2 - 8*r_um1 + 8*r_up1 - r_up2) ./ (12*hu);
    rv = (r_vm2 - 8*r_vm1 + 8*r_vp1 - r_vp2) ./ (12*hv);
    rq = (r_qm2 - 8*r_qm1 + 8*r_qp1 - r_qp2) ./ (12*hq);
    
    % Volume Jacobian
    % J = det([ru rv rq])
    %
    % Use absolute value because volume must be positive
    cross_uv = zeros(size(ru));

    cross_uv(:,1) = ru(:,2).*rv(:,3) - ru(:,3).*rv(:,2);
    cross_uv(:,2) = ru(:,3).*rv(:,1) - ru(:,1).*rv(:,3);
    cross_uv(:,3) = ru(:,1).*rv(:,2) - ru(:,2).*rv(:,1);

    detJ = cross_uv(:,1).*rq(:,1) + ...
           cross_uv(:,2).*rq(:,2) + ...
           cross_uv(:,3).*rq(:,3);

    jac = abs(detJ);

    w = wt .* jac;

    w = w(:);

end
