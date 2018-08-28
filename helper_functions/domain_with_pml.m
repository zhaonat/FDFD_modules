function [xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml)
%% Input parameters
% xrange: [xmin xmax], range of domain in x-direction without PML
% yrange: [ymin ymax], range of domain in y-direction without PML
% N: [Nx Ny], numbers of cells in x- and y-directions without PML
% Npml: [Nx_pml Ny_pml], numbers of cells in the x- and y-normal PML

%% Output parameters
% xrange: [xmin xmax], updated range of domain in x-direction including PML thickness
% yrange: [ymin ymax], updated range of domain in y-direction including PML thickness
% N: updated numbers of cells in the x- and y-directions including cells in PML
% dL: [dx dy], cell sizes
% Lpml: [Lx_pml Ly_pml], thicknesses of x- and y-normal PML

L = [diff(xrange) diff(yrange)];  % [Lx Ly]
dL = L./N;  % [dx dy]

Lpml = Npml .* dL;  % [Lx_pml, Ly_pml]
xrange = xrange + [-1 1] * Lpml(1);  % [xmin xmax] is updated
yrange = yrange + [-1 1] * Lpml(2);  % [ymin ymax] is updated

N = N + 2*Npml;  % [Nx Ny] is updated
