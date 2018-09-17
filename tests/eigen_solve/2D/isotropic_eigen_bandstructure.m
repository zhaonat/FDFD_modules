close all
clear
%% potential applications
% cavity mirrors: multi-frequency, lasers
% oled mirrors

%% design principle: we want magnitude of dielectric equal to magnitude of metallic dielectric?

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
wvlen = 2.2;  % wavelength in L0 
c0 = 3e8;
xrange = 0.1*[-1 1];  % x boundaries in L0
yrange = [-1 1];  % y boundaries in L0
a = diff(xrange);

thickness = 0.3;

L = [diff(xrange), diff(yrange)];
N = [150 210];  % [Nx Ny]
Npml = 1*[0 50];  % [Nx_pml Ny_pml]

Nx = N(1); Ny = N(2);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
omega_p = 3e15; gamma = 5.5e12;
omega = 2*pi*c0/wvlen*1e6;
epsilon_diel = 16;
epsilon_metal =  1; 1 - omega_p^2./(omega^2+1i*gamma*omega);
fill_factor = 0.2; %half metal, half dielectric

eps = hybrid_grating_grid(N, L, epsilon_diel, epsilon_metal,...
    fill_factor, thickness );

xlim = xrange;
ylim = [-thickness, thickness]/2;
figure()
visreal(eps, xrange, yrange);
colorbar;

neigs = 22;

%% MODE SCAN
band_structure = [];
wvlen_guess = 1;
figure()s;
plot(linspace(0, pi/L(1)*a/(2*pi), 200), linspace(0, pi/L(1)/(2*pi), 200))
hold on;
c = 1;
for Kx = linspace(0, pi/L(1), 20)
    %Kx =0.5*pi/L(1);
    K_vec=[Kx,0];
%     
%     [Hz, Ex, Ey, omega_eigs] = solveTE_BlochX(L0, wvlen_guess, xrange,...
%           yrange, eps, Kx, Npml, neigs);
    
    [Hz, Ex, Ey, omega_eigs,A] = eigensolve_TM(L0, wvlen_guess, xrange, ...
        yrange, eps, Npml, neigs, K_vec);
    %% we need to filter out spurious modes (particularly modes with 
    % heavily concentrated fields near/around the PML
    filtered_eigs = [];
%     for i = 1:n
%         scatter(Kx*a/(2*pi), (real(omega_eigs(i))*1e-6/(2*pi*c0)), '.b')
%         %text(Kx*a/(2*pi), (real(omega_eigs(i))*1e-6/(2*pi*c0)), num2str(i), 'Fontsize', 10);
%     end

    [filtered_modes, filtered_eigs] = ...
    mode_filtering(Hz, omega_eigs, eps, xlim, ylim, L, Npml);
    
    for i = 1:length(filtered_eigs)
        scatter(Kx*a/(2*pi), (real(filtered_eigs(i))*1e-6/(2*pi*c0)), '.b')
        %text(Kx*a/(2*pi), (real(omega_eigs(i))*1e-6/(2*pi*c0)), num2str(i), 'Fontsize', 10);
    end

    band_structure = [band_structure, omega_eigs];
    drawnow();
    hold on;
    c = c+1;

end

xlabel('Kx')
ylabel('eigen_frequency (s^-1)')



