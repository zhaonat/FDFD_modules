%% Dirichlet Solver
L0 = 1e-6;  % length unit: microns
wvlen = 1.55;  % wavelength in L0
diel = 1;
cellsize = 30;
inclusionSize = cellsize-10;

numcells = 1;
cellx = numcells; celly = numcells; cellz = 1;

Npml = [0 0 0];  % [Nx_pml Ny_pml]
xrange = numcells*[-1 1];  % x boundaries in L0
yrange = numcells*[-1 1];  % y boundaries in L0
zrange = numcells*[-1 1];
%[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  
% domain is expanded to include PML

%% Set up the permittivity.
SingleCellDims = cellsize*[1, 1, 1];
[eps_r, interiorCoords, borderCoords] = cubeDielectricGrid(cellx, celly,...
cellz, SingleCellDims,...
inclusionSize, diel);
N = size(eps_r);

%% Set up the current source density.
Mz = zeros(N); My = Mz; Mx = Mz;
ind_src = [1 1 1];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2), ind_src(3)) = 1;
JCurrentVector = [Mx; My; Mz];

%% implement an iterative solve for all of the Schur complement elements
%% Wonsoek's scalar parameter 1, -1, or 0
s = -1;

[A, b, Ao, bo, omega, c0, TepsSuper, TmuSuper] = ...
    solve3D_EigenEngine_Matrices_dirichlet(L0, wvlen, xrange, yrange, zrange, eps_r,...
    JCurrentVector, Npml, s);

%% solve
disp('iterative solution')
[x, flag, relres, iter, resvec] = qmr(A,b,1e-8, length(b));


%% fields
fieldsize = prod(N);
Ex = x(1:fieldsize);
Ex = reshape(Ex, N(1), N(2), N(3))

for i = 1:N(1)
   figure;
   visreal(1i*Ex(:,:,i), xrange, yrange) 
end

moviereal(Ex(:,:,5), xrange, yrange)



