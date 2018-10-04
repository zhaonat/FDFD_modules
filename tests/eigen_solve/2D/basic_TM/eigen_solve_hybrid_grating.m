close all
clear
%% potential applications
% this eigenstructure is NON-DISPERSIVE AND HENCE WRONG, but 
% it is still relatively illustrative

%% design principle: we want magnitude of dielectric equal to magnitude of metallic dielectric?

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
wvlen = 2.2;  % wavelength in L0 
c0 = 3e8;
xrange = 0.1*[-1 1];  % x boundaries in L0
yrange = [-1 1];  % y boundaries in L0
a = diff(xrange);


L = [diff(xrange), diff(yrange)]
N = [150 210];  % [Nx Ny]
Npml = 1*[0 40];  % [Nx_pml Ny_pml]

Nx = N(1); Ny = N(2);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
omega_p = 0.72*pi*1e15;%3e15;
gamma = 5.5e12;
omega = 2*pi*c0/wvlen*1e6;
epsilon_diel = 16;
epsilon_metal =  1 - omega_p^2./(omega^2+1i*gamma*omega);

fill_factor = 0.2; %half metal, half dielectric
thickness = 0.3;
eps = hybrid_grating_grid(N, L, epsilon_diel, epsilon_metal, fill_factor, thickness);
xlim = xrange;
ylim = [-thickness/2, thickness/2];
figure()
visabs(eps, xrange, yrange);
colorbar;

neigs = 51; %modes to get

%% MODE SCAN
band_structure = [];
wvlen_min = 1;
f = figure();
plot(linspace(0, pi/L(1)*a/(2*pi), 200), c0*1e6*linspace(0, pi/L(1)*1e-6/(2*pi*c0), 200))
hold on;
counter = 1;
for Kx = linspace(0, pi/L(1), 300)
    %Kx =0.5*pi/L(1);
    K_vec = [Kx,0];
    [Hz, Ex, Ey, omega_eigs,Am] = eigensolve_TM(L0, wvlen, xrange, ...
    yrange, eps, Npml, neigs, K_vec);
   
    %% we need to filter out spurious modes (particularly modes with 
    % heavily concentrated fields near/around the PML
    [filtered_modes, filtered_eigs, mask] = ...
    mode_filtering(Hz, omega_eigs,  eps, xlim, ylim, L, Npml);
    for j = 1:length(filtered_eigs)
        scatter(Kx*a/(2*pi),real(filtered_eigs(j))*1e-6/(2*pi*c0), '.b')
        %text(Kx*a/(2*pi), real(filtered_eigs(j))*1e-6/(2*pi*c0), num2str(real(filtered_eigs(j))*1e-6/(2*pi*c0)), 'Fontsize', 10);
    end
%     for j = 1:length(filtered_eigs)
%         visreal(filtered_modes{j}, xrange, yrange);
%         drawnow();
%     end    
    %plot(repmat(Kx*a/(2*pi), 1, length(filtered_eigs)), abs(filtered_eigs)*a*1e-6/(2*pi*c0), '.b');
    drawnow();
    hold on;
    counter = counter+1;
%     if(counter>1)
%         break
%     end

end

xlabel('Kx')
ylabel('eigen_frequency/(2*pi*c0)')
title('Kx band structure of the HG')
saveas(f, 'kx_band_structure_hybrid_grating.fig');



