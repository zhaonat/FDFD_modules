%% example of slab waveguide excitation
clear
close all

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 3.91;  % wavelength in L0

%% ==================== Create Data; Directory
Npml = [0,0];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
xrange = 2*[-0.5, 0.5];  % x boundaries in L0
yrange = 2*[-0.5, 0.5];  % y boundaries in L0

%% domain structuring
N = [100,100];  % [Nx Ny]
Nx = N(1); Ny = N(2);
N0 = N;
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml); 

h = dL(1);
%% Set up the magnetic current source density.
xa = 1:N(1);
ya = 1:N(2);
[Y,X] = meshgrid(ya,xa);

k0 = (2*pi)/wvlen;
x0 = [0.1, 0.1]; %% initial guess
omega = k0;
d = 0.5;
eps_diel = 12;
eps_air = 1;
[k_soln,fval] = fsolve(@(x) dielectricSlabEq(x, omega,d, eps_diel, eps_air), x0)

kappa = k_soln(1);
k = k_soln(2);
right = d; left = -d;
slab_range = [-d, d];
Jz_line = slab_mode_excitation(xrange, h, k, kappa, right, left);
figure;
plot(Jz_line)
Mz = zeros(N);
Mz(:,90) = Jz_line;

%% slabe waveguide structure
eps = ones(N);
eps(25:75,:) = 12;

% %% coupling
% eps(75:95, :) = 12;

%% SOLVERS
[Hz, Ex, Ey, A, omega_sim,b] = solveTM(L0,wvlen, xrange, yrange, eps, Mz, Npml);

figure; 
visreal(1i*Hz, xrange, yrange)

figure; imagesc(eps)

    

    

