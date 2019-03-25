%% test script for TM Ex Ey eigensolve

close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
mu0 = 4*pi*1e-7*L0;
eps0 = 8.854e-12*L0;
c0 = 1/sqrt(mu0*eps0);
xrange = 0.5*[-1,1];  % x boundaries in L0
yrange = 0.5*[-1,1];  % y boundaries in L0
a = diff(xrange);
L = [diff(xrange), diff(yrange)];
N = [120 120];  % [Nx Ny]
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

Ky = 0;
neigs = 4;
kx_guess = 0;
%wvlen_scan = linspace(2, 20, 400);
wvlen_scan = logspace(log10(0.9), log10(10), 10000);
kx_spectra = [];
% c = 1;
% for wvlen = (wvlen_scan)
%     omega = 2*pi*c0/wvlen;
%     [Ez_modes, Hx_modes, Hy_modes, eigenvals] = ...
%         eigensolve_TE_dispersive_Kx(L0, omega, xrange, ...
%         yrange, eps_r, Npml, neigs, kx_guess, Ky);
%     kx_spectra = [kx_spectra, eigenvals];
% %     if(mod(c,20) == 0)
% %        figure()
% %        visreal(Ez_modes{1}, xrange, yrange)
% %        title(c)
% %        drawnow();
% %     end
%     c = c+1;    
% end


parfor c = 1:length(wvlen_scan)
    wvlen = wvlen_scan(c);
    omega = 2*pi*c0/wvlen;
    [Ez_modes, Hx_modes, Hy_modes, eigenvals] = ...
        eigensolve_TE_dispersive_Kx(L0, omega, xrange, ...
        yrange, eps_r, Npml, neigs, kx_guess, Ky);
    kx_spectra(c,:) =  eigenvals;
end


omega_scan = (2*pi*c0./wvlen_scan)*(a/(2*pi*c0));
figure();
plot(real(kx_spectra)*a, omega_scan, '.b');
hold on;
plot(imag(kx_spectra), omega_scan, '.r');
ylabel('omega*a/2 pi c0')
xlabel('ka')
xlim([-pi, pi])

figure();
plot(real(kx_spectra)*a, omega_scan, '.b');
hold on;
ylabel('omega*a/2 pi c0')
xlabel('ka')
xlim([-pi, pi])
ylim([0,0.8])
save('photonic_circle_bandstructure_for_comparison_with_Johannopoulos_book.mat')


