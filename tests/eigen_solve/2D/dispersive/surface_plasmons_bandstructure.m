%% test script for TM Ex Ey eigensolve

close all
clear

% for surface plasmons, you would want a PEC (no propagating waves from
% spps)

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
mu0 = 4*pi*1e-7*L0;
eps0 = 8.854e-12*L0;
c0 = 3e8;
xrange = 0.2*[-1,1];  % x boundaries in L0
yrange = 1*[-1,1];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [100 150];  % [Nx Ny]
Npml = 1*[0 10];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);
wvlen = 2;
%% Set up the permittivity.
omega = 2*pi*c0/(wvlen*L0);
omega_p =  0.72*pi*1e15;
gamma = 5.5e12;
epsilon_metal = 1-omega_p^2/(omega^2-1i*gamma*omega);
epsilon = ones(N);
x = 1:N(1);
y = 1:N(2);
[xx, yy] = meshgrid(x,y);
half_ny = 14;
xbounds = [-0.1, 0.1];
ybounds = [-half_ny*dL(2), half_ny*dL(2)];
epsilon(:, 1:cy) = epsilon_metal;
figure();
visreal(epsilon, xrange, yrange);
drawnow();

%% eigensolve
neigs = 100;
kx_guess = 0*pi/L(1);
[Hz_modes, Ex_modes, Ey_modes, kx_eigs] = ...
    eigensolve_TM_dispersive_Kx(L0, omega, xrange, yrange, epsilon, Npml, neigs, kx_guess);
[filtered_modes, filtered_k] = ...
    mode_filtering(Hz_modes, kx_eigs, epsilon, xbounds, ybounds, L, Npml);
filtered_modes = Hz_modes; filtered_k = kx_eigs;
for i = 1:length(filtered_k)
    figure();
    Kx = filtered_k(i);
    visreal(filtered_modes{i}, xrange, yrange);
    title(strcat(num2str(i), ', ', num2str(real(Kx)/(2*pi)*diff(xrange))));  
end
