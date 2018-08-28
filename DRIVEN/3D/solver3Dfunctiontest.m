%% Solver3D function test

c0 = 3*10^8;
L0 = 1e-6;  % length unit: microns
wvlen = 4.0*L0;  % wavelength in L0
xrange = [-5 5]*L0;  % x boundaries in L0
yrange = [-5 5]*L0;  % y boundaries in L0
zrange = [-5 5]*L0;
N = [20 20 20];  % [Nx Ny]
Npml = [0 0 0];  % [Nx_pml Ny_pml]
mu0 = 4*pi*10^-7; mu_0 = mu0; mu = mu0;
eps0 = 8.85*10^-12; eps_0 = eps0;
%[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
eps_r = ones(N);

%% Set up the 1agnetic current source density.
Mz = zeros(N); My = Mz; Mx = Mz;
ind_src = [10 10 10];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2), ind_src(3)) = 1;
JCurrentVector = [Mx; My; Mz];

%% Execute Function Solver
[Hz, Ex, Ey, A] = solve3D(wvlen, xrange, yrange, zrange, eps_r, JCurrentVector, Npml);
