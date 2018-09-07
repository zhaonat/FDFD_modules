clear all; close all; clc;

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
wvlen = 5.0;  % wavelength in L0
xCells = 40, yCells = 30;
xrange = xCells*[-5 5];  % x boundaries in L0
yrange = yCells*[-5 5];  % y boundaries in L0
Npml = [0 0];  % [Nx_pml Ny_pml]

%% Set up the permittivity.
eps_air = multiRandomCellDielectric(xCells, yCells, 100, 100, Npml);
N = size(eps_air);  % [Nx Ny]
Nx = N(1) 
Ny = N(2)
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;

%% Solve TE equations.
tic
[Hz, Ex, Ey,A,b] = solveTE(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);
toc

%% Visualize the solution.
figure;
visabs(Hz, xrange, yrange)

%%
figure;
visabs(Ex, xrange, yrange - 0.5*dL(2));

figure;
visabs(Ey, xrange - 0.5*dL(1), yrange);

%% Show the movie of the oscillating field.
% figure;
% moviereal(Hz, xrange, yrange)
