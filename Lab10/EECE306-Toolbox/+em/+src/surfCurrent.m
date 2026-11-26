function s = surfCurrent(K,surf,uspan,vspan,Nu,Nv,rule)
%SURFCURRENT Discretize a surface-current density.
%   K is a function handle returning one 1-by-3 current-density vector
%   (A/m) at each surface position. Element vectors are K(r_i)*dS_i.
if nargin < 7 || isempty(rule), rule = 'midpoint'; end
[pos,w] = em.quad.surf(surf,uspan,vspan,Nu,Nv,rule);
if ~isa(K,'function_handle')
    error('em.src.surfCurrent:InvalidK', ...
        'K must be a function handle returning 1-by-3 vectors.');
end
try
    kval = K(pos);
catch
    kval = [];
end
if numel(kval) ~= 3*size(pos,1)
    kval = zeros(size(pos));
    for k = 1:size(pos,1)
        kk = K(pos(k,:));
        if numel(kk) ~= 3
            error('K must return one 1-by-3 vector per position.');
        end
        kval(k,:) = reshape(kk,1,3);
    end
else
    kval = reshape(kval,[],3);
end
if ~isequal(size(kval),size(pos)) || any(~isfinite(kval(:)))
    error('K must return finite 1-by-3 vectors at every position.');
end
s.type = 'current';
s.pos = pos;
s.w = w(:);
s.q = [];
s.Idl = kval .* w(:);
s.label = '';
end
