close all
clear
%% Dirichlet Solver
L0 = 1e-6;  % length unit: microns
wvlen = 1.55;  % wavelength in L0
diel = 1;
cellsize = 30;
inclusionSize = cellsize-10;
N = [20,20,20]

numcells = 1;
cellx = numcells; celly = numcells; cellz = 1;

Npml = 1*[3,3,3];  % [Nx_pml Ny_pml]
xrange = numcells*[-1 1];  % x boundaries in L0
yrange = numcells*[-1 1];  % y boundaries in L0
zrange = numcells*[-1 1];
N = N+Npml;
%[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  
% domain is expanded to include PML
%% Set up the permittivity.
SingleCellDims = cellsize*[1, 1, 1];
eps_r = ones(N);
N = size(eps_r);
dL = 2./N;

%% Set up the current source density.
Mz = zeros(N); My = Mz; Mx = Mz;
ind_src = [1 1 1];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2), ind_src(3)) = 1;
JCurrentVector = [Mx; My; Mz];

%% implement an iterative solve for all of the Schur complement elements
%% Wonsoek's scalar parameter 1, -1, or 0
s = -1;

[A, b, Ao, bo, omega, c0, TepsSuper, TmuSuper, ...
    Sxf, Syf, Szf, Sxb, Syb, Szb] = ...
    solve3D_EigenEngine_Matrices_dirichlet(L0, wvlen, xrange, yrange, zrange, eps_r,...
    JCurrentVector, Npml, s);

%% SYMMETRIZER
dx = dL(1); dy = dL(2); dz = dL(3);
%% stretch factor
sxf = diag(Sxf);
syf = diag(Syf);
szf = diag(Szf);

numerator = sqrt((sxf.*syf.*szf));
numerator = vertcat(numerator, numerator, numerator);
denominator = 1./(numerator);
dim = prod(N)*3
Pr = spdiags(numerator, 0, dim, dim);
Pl = spdiags(denominator, 0, dim, dim);

SA = Pl^-1*(A*Pr^-1);

%% test of the QMR
disp('solve')
tic
[x, flag, relres, iter, resvec] = qmr(A,b, 1e-6, length(b));
toc
tic
[x2, flag2, relres2, iter2, resvec2] = qmr(SA,Pl\b, 1e-6, length(b));
toc
x2 = Pr\x2;

