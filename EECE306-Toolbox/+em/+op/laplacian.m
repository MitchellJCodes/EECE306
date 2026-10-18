function L = laplacian(f, r, h)
%LAPLACIAN Pointwise Laplacian of a Scalar Field
%   L = EM.OP.LAPLACIAN(F,R,H) evaluates d2f/dx2+d2f/dy2 at rows of r.
%   f accepts either an N-by-3 point array (like em.op.grad) or separate
%   x,y arrays. r is N-by-2 or N-by-3. h is a positive scalar or [hx hy].
%   The second derivative uses centred second differences.
if nargin < 3 || isempty(h), h = 1e-4; end
if ~isa(f,'function_handle'), error('f must be a function handle'); end
if size(r,2) ~= 2 && size(r,2) ~= 3
    error('r must be an N-by-2 or N-by-3 coordinate array');
end
if isscalar(h)
    hx = h; hy = h;
elseif isnumeric(h) && numel(h) == 2
    hx = h(1); hy = h(2);
else
    error('h must be a positive scalar or [hx hy]');
end
if any(~isfinite([hx hy])) || any([hx hy] <= 0)
    error('h must contain positive finite step sizes');
end

rp = r; rm = r;
rp(:,1)=rp(:,1)+hx; rm(:,1)=rm(:,1)-hx;
fp = evalField(f,rp); fm = evalField(f,rm); f0 = evalField(f,r);
d2x = (fp - 2*f0 + fm)/(hx^2);
rp = r; rm = r;
rp(:,2)=rp(:,2)+hy; rm(:,2)=rm(:,2)-hy;
fp = evalField(f,rp); fm = evalField(f,rm);
d2y = (fp - 2*f0 + fm)/(hy^2);
L = d2x+d2y;
end

function y = evalField(f,r)
try
    y = f(r);
catch
    y = f(r(:,1),r(:,2));
end
if isscalar(y) && size(r,1) > 1, y = repmat(y,size(r,1),1); end
y = y(:);
if numel(y) ~= size(r,1) || any(~isfinite(y))
    error('f must return one finite scalar per observation point');
end
end

