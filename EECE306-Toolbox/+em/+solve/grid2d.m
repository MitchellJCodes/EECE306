function G = grid2d(xlim, ylim, Nx, Ny)
%GRID2D 2D Grid Struct
%   G = GRID2D(XLIM, YLIM, NX, NY) returns a returns a struct with G.x,
%   G.y, G.X, G.Y, G.hx and G.hy.

x = linspace(xlim(1), xlim(2), Nx);
y = linspace(ylim(1), ylim(2), Ny);
[X, Y] = meshgrid(x, y);          % Ny-by-Nx
G = struct('x', x, 'y', y, 'X', X, 'Y', Y, ...
           'hx', x(2)-x(1), 'hy', y(2)-y(1), 'Nx', Nx, 'Ny', Ny);
end
