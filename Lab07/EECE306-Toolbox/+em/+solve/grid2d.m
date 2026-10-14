function G = grid2d(xlim, ylim, Nx, Ny)
%GRID2D Construct a uniform Cartesian grid for a 2-D boundary-value problem.
%   G = em.solve.grid2d(xlim,ylim,Nx,Ny) creates Nx points in x and Ny
%   points in y, including both endpoints. G.X and G.Y are Ny-by-Nx arrays
%   in meshgrid layout; hx and hy are the uniform spacings.

if nargin ~= 4
    error('grid2d requires xlim, ylim, Nx, and Ny');
end
if numel(xlim) ~= 2 || numel(ylim) ~= 2 || ...
        any(~isfinite([xlim(:); ylim(:)])) || ...
        xlim(2) <= xlim(1) || ylim(2) <= ylim(1)
    error('xlim and ylim must be finite increasing two-element intervals');
end
if ~isscalar(Nx) || ~isscalar(Ny) || Nx < 3 || Ny < 3 || ...
        Nx ~= floor(Nx) || Ny ~= floor(Ny)
    error('Nx and Ny must be integers greater than or equal to 3');
end

G.x = linspace(xlim(1), xlim(2), Nx);
G.y = linspace(ylim(1), ylim(2), Ny);
[G.X, G.Y] = meshgrid(G.x, G.y);
G.hx = G.x(2) - G.x(1);
G.hy = G.y(2) - G.y(1);
G.Nx = Nx;
G.Ny = Ny;
G.xlim = xlim(:).';
G.ylim = ylim(:).';
end
