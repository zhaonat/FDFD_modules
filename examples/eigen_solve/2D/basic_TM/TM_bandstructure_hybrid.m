close all
clear
%% potential applications
% cavity mirrors: multi-frequency, lasers
% oled mirrors

%% design principle: we want magnitude of dielectric equal to magnitude of metallic dielectric?

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;
xrange = 0.1*[-1 1];  % x boundaries in L0
yrange = [-1 1];  % y boundaries in L0
a = diff(xrange);

thickness = 0.3;

L = [diff(xrange), diff(yrange)];
N = [100 120];  % [Nx Ny]
Npml = 1*[0 10];  % [Nx_pml Ny_pml]

Nx = N(1); Ny = N(2);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML



xlim = xrange;
ylim = [-thickness, thickness]/2;

neigs = 100;

%% MODE SCAN
band_structure = [];
wvlen_guess = 1;
figure();
plot(linspace(0, pi/L(1)*a/(2*pi), 200), linspace(0, pi/L(1)/(2*pi), 200))
hold on;
c = 1;
omega_p = 3e15; gamma = 5.5e12;

for omega = linspace(omega_p/5, 1.2*omega_p, 40)
    wvlen = 2*pi*c0/(omega)/1e-6;
    kx_guess = -0.1*pi/wvlen;

    %% Set up the permittivity.
    epsilon_diel = 16;
    epsilon_metal =  1 - omega_p^2./(omega^2-1i*gamma*omega);
    fill_factor = 0.2; %half metal, half dielectric

    eps = hybrid_grating_grid(N, L, epsilon_diel, epsilon_metal,...
        fill_factor, thickness );
    
    [Hz_modes, Ex_modes, Ey_modes, kx_eigs] = ...
        eigensolve_TM_dispersive_Kx(L0, omega, xrange, yrange,...
        eps, Npml, neigs, kx_guess);
    [filtered_modes, filtered_k] = ...
        mode_filtering(Hz_modes, kx_eigs, eps, xlim, ylim, L, Npml);

    for i = 1:length(filtered_k)
        scatter(abs(real(filtered_k(i))*a/(2*pi)), omega*1e-6/(2*pi*c0), '.b')
        %text(Kx*a/(2*pi), (real(omega_eigs(i))*1e-6/(2*pi*c0)), num2str(i), 'Fontsize', 10);
    end

    band_structure = [band_structure, kx_eigs];
    drawnow();
    hold on;
    c = c+1;

end

xlabel('Kx')
ylabel('eigen_frequency (s^-1)')



