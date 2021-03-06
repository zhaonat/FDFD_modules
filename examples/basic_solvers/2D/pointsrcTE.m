%clear all; close all; clc;

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
wvlen = 4;  % wavelength in L0
xrange = [-5 5];  % x boundaries in L0
yrange = [-5 5];  % y boundaries in L0
N = [300 220];  % [Nx Ny]
Npml = 0*[15 15];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0

%% Note on grid resolution of the system
% dx/(wvlen) ~1/20 or smaller

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
%FINAL GRID PARAMETERS ARE DETERMINED AT THIS POINT
%xrange and yrange are slightly larger


%% NOTE dL is not in SI units when it comes out of domain_with_pml; for our purposes, that is okay
%for now
resolutionFactor = max([dL(1)/N(1) dL(2)/N(2)]); %dx/N ~meters
%spatially, what is the smallestlength scale that we have to resolve

%% Set up the permittivity.
eps_air = ones(N);


%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;

%% Solve TE equations.
[Hz, Ex, Ey,A] =...
    solveTE(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);

%% Visualize the H-field solution.
figure;
visabs(Hz, xrange, yrange)

%% visualize E fields
figure;
visabs(Ex, xrange, yrange - 0.5*dL(2));

figure;
visabs(Ey, xrange - 0.5*dL(1), yrange);

%% Show the movie of the oscillating field.
figure;
moviereal(Hz, xrange, yrange);


%% Visualize the H-field solution.
figure;
visreal(1j*Hz, xrange, yrange)

%% visualize E fields
figure;
visreal(Ex, xrange, yrange - 0.5*dL(2));

figure;
visreal(Ey, xrange - 0.5*dL(1), yrange);


