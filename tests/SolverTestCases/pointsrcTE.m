clear all; close all; clc;

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
wvlen = 8;  % wavelength in L0
xrange = [-5 5];  % x boundaries in L0
yrange = [-5 5];  % y boundaries in L0
N = [100 100];  % [Nx Ny]
Npml = 0*[40 40];  % [Nx_pml Ny_pml]
Nx = N(1); Ny = N(2);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
epsilon = 1;
numCells = 10;
cellsize = Nx/numCells;
eps_air = ones(N);
eps_air(30:70, 30:70) = 12;
%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;

%% Solve TE equations.
tic
[Hz, Ex, Ey,A,omega,b] = solveTM(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);
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
