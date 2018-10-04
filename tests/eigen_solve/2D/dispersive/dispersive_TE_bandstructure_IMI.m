%% test script for TM Ex Ey eigensolve

close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
mu0 = 4*pi*1e-7*L0;
eps0 = 8.854e-12*L0;
c0 = 1/sqrt(eps0*mu0);
xrange = [-1,1];  % x boundaries in L0
yrange = [-2,2];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [80 80];  % [Nx Ny]
Npml = 1*[0 0];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);
wvlen = 3;
wvlen_scan = logspace(log10(2),log10(10),30);
figure();
omega_p =  0.72*pi*1e15;
gamma = 5.5e12;

for wvlen = wvlen_scan
    %% Set up the permittivity.
    omega = 2*pi*c0/(wvlen*L0);

    epsilon_metal = 16; %1-omega_p^2/(omega^2-1i*gamma*omega);
    epsilon_diel = 16; %epsilon_metal;
    fill_factor = 0.2;
    thickness = 0.6;
    eps = ones(N);
    num_cells = 1;
    y_grid_center = L(2)/2;
    y_center = y_grid_center-0.3;
    eps = hybrid_grating_multi_unit_cell_add(eps,num_cells, N, L, epsilon_diel,...
        epsilon_metal, fill_factor, thickness, y_center);
    y_center = y_grid_center+0.3;
    eps = hybrid_grating_multi_unit_cell_add(eps,num_cells, N, L, epsilon_diel,...
        epsilon_metal, fill_factor, thickness, y_center);
    xlim = xrange;
    ylim = [-1, 1];

    %% eigensolve
    neigs = 30;
    kx_guess = 0.25*pi/L(1);
    [Ez_modes, Hx_modes, Hy_modes, kx_eigs] = ...
        eigensolve_TE_dispersive_Kx(L0, omega, xrange, yrange, eps, Npml, neigs, kx_guess);
    [filtered_modes, filtered_k, mask] = ...
        mode_filtering(Ez_modes, kx_eigs, eps, xlim, ylim, L, Npml);

    for i = 1:length(filtered_k)
        scatter( abs(filtered_k(i))*L(1)/(2*pi), omega, '.b');
        hold on;
 
    end
    drawnow();
end
kx_scan = linspace(0, pi/L(1), 1000);
omega_light = c0*kx_scan;
hold on;
plot(kx_scan, omega_light);

for i = 1:length(filtered_k)
    figure();
    Kx = filtered_k(i);
    visreal(filtered_modes{i}, xrange, yrange);
    title(strcat(num2str(i), ', ', num2str(real(Kx)/(2*pi)*diff(xrange))));  
end



% figure();
% moviereal(filtered_modes{2}, xrange, yrange)

