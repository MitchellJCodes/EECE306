function C = capacitance(G,epsr,mask,Vfix,opts)
%CAPACITANCE Compute capacitance per unit length from stored electric energy.
% C = 2W/V0^2, where V0 is the difference between fixed conductor potentials.
if nargin<5 || isempty(opts), opts=struct(); end
[V,~]=em.solve.poisson2d(G,0,epsr,mask,Vfix,opts);
if isscalar(mask), mask=false(size(G.X)); end
if isscalar(Vfix), Vfix=repmat(Vfix,size(G.X)); end
vals=Vfix(logical(mask));
if numel(vals)<2, error('capacitance:conductors','mask must include at least two fixed-potential nodes.'); end
V0=max(vals)-min(vals);
if ~isfinite(V0) || V0<=0, error('capacitance:voltage','Fixed conductors must have different potentials.'); end
W=em.derive.energyE(G,V,epsr,opts);
C=2*W/(V0^2);
end
