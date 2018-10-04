%close all
%clear
%% potential applications
% cavity mirrors: multi-frequency, lasers

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;

xrange = [-2 2];  % x boundaries in L0
yrange = 1*[-2 2];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [800 800];  % [Nx Ny]
Npml = 1*[40 40];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] =domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);

%% Set up the permittivity
wvlen_scan = linspace(1,2.6, 20);

%parfor wvlen = wvlen_scan
reflected_power = [];
nx = 0;
%wvlen = 2.2632;
wvlen = 1.7;
k0 = 2*pi/wvlen;
omega_p = 0.72*pi*1e15;%3e15; %omega_p was 3 which gave us good results...
gamma = 400e12; %20e12; % (5.5e12 is the default)
omega = 2*pi*c0/wvlen*1e6;
epsilon_diel = 16;
epsilon_metal =  1 - omega_p^2./(omega^2-1i*gamma*omega); %MINUS SIGN FOR FDFD
epsilon_diel = epsilon_metal;
thickness = 0.15;
fill_factor = 0.1; %half metal, half dielectric

delta_arc = 6*pi/180;
inner_radius = 0.5; outer_radius = 0.7;
eps = ones(N);
eps = curved_stripe(eps, N,xrange, yrange, ...
    inner_radius, outer_radius, delta_arc, epsilon_metal, epsilon_diel);

delta_arc_2 = 3*pi/180;
inner_rad_2 = 1.3; outer_rad_2 = 1.5;
eps = curved_stripe(eps, N,xrange, yrange, ...
    inner_rad_2, outer_rad_2, delta_arc_2, epsilon_metal, epsilon_diel);

figure(); visabs(eps, xrange, yrange);
drawnow();

%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = [floor(N(2)/2), floor(N(2)/2)];  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(250, 450) = 1;

%% set up electric current source density
Jx = zeros(N);
Jy = zeros(N);
% Jx(ind_src(1), ind_src(2)) = 1;
% Jx(ind_src(1), ind_src(2)+1) = -1;

%% Solve TE (photonic), TM (practical), depends on convention equations.
tic
[Hz, Ex, Ey] = ...
    solveTE(L0, wvlen, xrange, yrange, eps, Mz, Npml);
toc

visabs(Hz, xrange, yrange)
figure();
visabs(Hz/max(max(abs(Hz))), xrange, yrange)



