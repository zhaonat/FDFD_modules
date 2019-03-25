%% test script for TM Ex Ey eigensolve
clear all
close all

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
mu0 = 4*pi*1e-7*L0;
eps0 = 8.854e-12*L0;
c0 = 1/sqrt(mu0*eps0);
xrange = 0.5*[-1,1];  % x boundaries in L0
yrange = 0.5*[-1,1];  % y boundaries in L0
a = diff(xrange);
L = [diff(xrange), diff(yrange)];
N = [150 150];  % [Nx Ny]
Npml = 1*[0 0];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);
i = 1:Nx;
j = 1:Ny;

[I,J] = meshgrid(i,j);
eps_r = ones(Nx,Ny);
radius = 0.2;
R = (radius/a)*Nx;

eps_r((J-cx).^2+(I-cy).^2 <= R^2) = 8.9;

visreal(eps_r, xrange, yrange)
drawnow();

neigs = 30;
kx_guess = 0;
%wvlen_scan = linspace(2, 20, 400);
kx_spectra = [];

k_scan = linspace(-pi/a, pi/a, 240);
parfor c = 1:length(k_scan)
    kx = k_scan(c);
    wvlen = abs(2*pi/kx);
%     omega = 2*pi*c0/wvlen;
    K_vec = [kx,0];
    [Ez_modes, Hx_modes, Hy_modes, eigenvals] = eigensolve_TE(L0, wvlen, xrange, ...
    yrange, eps_r, Npml, neigs, K_vec);
    kx_spectra(c,:) =  eigenvals;
end

%continuation scan
parfor c = 1:length(k_scan)
    kx = 0+1i*k_scan(c);
    wvlen = abs(2*pi/kx);
%     omega = 2*pi*c0/wvlen;
    K_vec = [kx,0];
    [Ez_modes, Hx_modes, Hy_modes, eigenvals] = eigensolve_TE(L0, wvlen, xrange, ...
    yrange, eps_r, Npml, neigs, K_vec);
    cont_spectra(c,:) =  eigenvals;
end

figure();
plot(k_scan, real(kx_spectra)*a/(2*pi*c0), '.b');
hold on;
plot(k_scan, real(cont_spectra)*a/(2*pi*c0), '.r');
hold on;
ylabel('omega*a/2 pi c0')
%xlim([-0.5*pi/a, 0.5*pi/a])
xlabel('ka')

save('non_dispersive_photonic_circle_bandstructure_for_comparison_with_Johannopoulos_book.mat')


