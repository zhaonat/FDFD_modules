%% vacuum dipole test
close all
clear
L0 = 1e-6;  % length unit: microns
wvlen = 2;  % wavelength in L0
diel = 1;
%max cell size can easily be 100
cellsize = 50;

numcells = 1;
cellx = numcells; celly = numcells; cellz = 1;
N = [cellsize, cellsize, cellsize];
Npml = 2*[5,5,5];  % [Nx_pml Ny_pml]
xrange = numcells*[-1 1];  % x boundaries in L0
yrange = numcells*[-1 1];  % y boundaries in L0
zrange = numcells*[-1 1];

%[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  

%% Set up the permittivity.
eps_r = ones(N);

%% Set up the current source density.
Mz = zeros(N); My = Mz; Mx = Mz;
ind_src = ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2), ind_src(3)) = 1;
JCurrentVector = [Mx; My; Mz];

%% implement an iterative solve for all of the Schur complement elements
%% Wonsoek's scalar parameter 1, -1, or 0, -1 is the best choice
s = -1;
disp('build system')
[A, b, Ao, bo, omega, TepsSuper, TmuSuper] = ...
    solve3D_matrices(L0, wvlen, xrange, yrange, zrange, eps_r,...
    JCurrentVector, Npml, s);

%% solve
disp('iterative solution')
[x, flag, relres, iter, resvec] = qmr(A,b,1e-8, length(b));

%% fields
fieldsize = prod(N);
Ex = x(1:fieldsize);
Ex = reshape(Ex, N(1), N(2), N(3));

for i = 1:N(1)
   figure;
   visabs(Ex(:,:,i), xrange, yrange) 
   caxis([0,0.02])
end
figure();
moviereal(Ex(:,:,35), xrange, yrange)



