function B_r = B(s,r,mur)
%B Magnetic flux density from a current source

if nargin < 3 || isempty(mur), mur = 1; end
if ~isscalar(mur) || ~isfinite(mur) || mur < 0
    error('mur must be a finite nonnegative scalar.');
end
B_r = em.const.mu0() .* mur .* em.field.H(s,r);
end
