%% example of slab waveguide excitation
clear
close all

%% ================ ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 3.91;  % wavelength in L0
mu0 = 4*pi*1e-7*L0;
eps0= 8.854e-12*L0;
c0 = 1/sqrt(mu0*eps0);

%% ==================== Create Data; Directory
Npml = [0,10];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0
xrange = 4*[-0.5, 0.5];  % x boundaries in L0
yrange = 6*[-0.5, 0.5];  % y boundaries in L0

%% domain structuring
N = [200,600];  % [Nx Ny]
Nx = N(1); Ny = N(2);
N0 = N;
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml); 

h = dL(1);
%% Set up the magnetic current source density.
xa = 1:N(1);
ya = 1:N(2);
[Y,X] = meshgrid(ya,xa);

k0 = (2*pi)/wvlen;
kz_guess = k0;
omega = c0*k0;

eps_diel = 12; e_slab = eps_diel;
eps_air = 1; e_out = eps_air;
mode = 'TE';
d = 0.5; %d = 1/2 of the actual slab thickness
kz_max = sqrt(omega^2*(eps0*mu0)*(e_slab-e_out));

fun = @(kz)tan(kz*d) - (sqrt((kz_max*d).^2-(kz*d).^2))./(kz*d);

%plot this equation;
kz_scan = linspace(0,3*k0,100);
figure(); plot(kz_scan, fun(kz_scan));

[kz,fval] = fzero(fun, kz_guess);

fun(kz)
kappa = sqrt(kz_max^2-kz^2);
k = kz;
right = d; left = -d;
slab_range = [-d, d];

%% requires us to have solved k and kappa
Jz_line = slab_mode_excitation(xrange, h, k, kappa, right, left);
figure;
x = linspace(xrange(1), xrange(2), N(1));
plot(x,Jz_line)
Mz = zeros(N);
Mz(:,90) = Jz_line;
figure();
visreal(Mz, xrange, yrange);
drawnow();

%% slabe waveguide structure
eps = ones(N);
eps(75:125,:) = 12;


%% SOLVERS
[Ez, Hx, Hy, A, omega_sim,b] = solveTM(L0,wvlen, xrange, yrange, eps, Mz, Npml);

figure; 
visreal(1i*Ez, xrange, yrange)

figure; imagesc(eps)

    

    

