%ptsrc_TE with a dirichlet boundary condition encoded into DWS
% clear all; close all; clc;

%% 10/31/2018 Currently, we know that the main issue is the on-diagonal, which we cannot match
%% to the PEC_PMC MASK RESULT, the diagonals of the dirichlet solution consistently underestimate
%% the mask solution...

%% Parameter Setup
L0 = 1e-6;  % length unit: microns
wvlen = 2;  % wavelength in L0
xrange = 1.5*[-1 1];  % x boundaries in L0
yrange = 1.5*[-1 1];  % y boundaries in L0
N = 30;
Npml = 0*[5 5];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
eps_air = ones(N);

%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;

%% Solve TE equations.
[Hz, Ex, Ey, A, omega,b, Sxf, Dxf, Dxb, Dyf, Dyb] =...
    solveTE_dirichlet(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);

%% Visualize the solution.
figure;
visabs(Hz, xrange, yrange)

%%
figure;
visabs(Ex, xrange, yrange - 0.5*dL(2));

figure;
visabs(Ey, xrange - 0.5*dL(1), yrange);

%% Show the movie of the oscillating field.
figure;
moviereal(Hz, xrange, yrange)




