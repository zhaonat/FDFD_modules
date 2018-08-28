clear all; close all; clc;

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
wvlen = 1;  % wavelength in L0
numCells = 1;
xrange = numCells*[-1 1];  % x boundaries in L0
yrange = numCells*[-1 1];  % y boundaries in L0
N = [80 80];  % [Nx Ny]
Npml = 0*[15 15];  % [Nx_pml Ny_pml]
Nx = N(1); Ny = N(2);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
cellSize = Nx/numCells; epsilon = 1;
eps_air = multiRandomCellDielectric(numCells,numCells, cellSize, cellSize, Npml, epsilon);
%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;

%% Solve TE equations.
tic
[Hz, Ex, Ey,A,b] = solveTE_dirichlet(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);
toc
tic
[Hz2, Ex2, Ey2,A2,b2] = solveTE(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);
toc

%% Visualize the solution.
figure;
visreal(Hz, xrange, yrange)

%%
figure;
visreal(Ex, xrange, yrange - 0.5*dL(2));

figure;
visreal(Ey, xrange - 0.5*dL(1), yrange);

%% Show the movie of the oscillating field.
figure;
moviereal(Hz, xrange, yrange)
