clear all; close all; clc;

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
wvlen = 3;  % wavelength in L0
xrange = [-5 5];  % x boundaries in L0
yrange = [-5 5];  % y boundaries in L0
N = [80 80];  % [Nx Ny]
Npml = [10 10];  % [Nx_pml Ny_pml]
Nx = N(1); Ny = N(2);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
epsilon = 12;
numCells = 1;
cellsize = Nx/numCells;
featureDims = [cellsize/4, cellsize/8, cellsize/4];
eps_air = multiRandomCellDielectric(numCells, numCells, cellsize, ...
    cellsize, Npml, epsilon, featureDims);
%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;

%% Solve TE equations.
tic
[Hz, Ex, Ey,A,omega,b,  Sxf, Syf, Sxb, Syb] = solveTE(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);
toc

%% SOLVER CODE HERE
N = size(Sxf);

dx = dL(1); dy = dL(2);
sxf = diag(Sxf);
syf = diag(Syf);
% sxf = 1-imag(diag(Sxf));
% syf = 1-imag(diag(Sxb));
% sxf = (diag(Sxf));
% syf = (diag(Sxb));
% sxf(sxf>1) = sxf(sxf>1)*1i;
% syf(syf>1) = syf(syf>1)*1i;

numerator = sqrt(sxf).*sqrt(syf);
denominator = 1./(numerator);
Pr = spdiags(numerator, 0, N(1), N(2));
Pl = spdiags(denominator, 0, N(1), N(2));

%[Pl, Pr] = SCSymmetrizer2D(Sxf, Syf, Sxb, Syb, dL);
sym = Pl^-1*A*Pr^-1;

sym(1:5, 1:5)
A(1:5, 1:5)

tic
x = (sym\(Pl\b));
x = Pr\x;
Hz2 = reshape(x,100, 100);
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
