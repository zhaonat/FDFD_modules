

clear
close all

c0 = 3*10^8;
L0 = 1e-6;  % length unit: microns
R = [];
for wvlen = 2:0.1:3
    xrange = 3*[-0.25, 0.25];  % x boundaries in L0
    yrange = [-5 5];  % y boundaries in L0
    N = [500 300];  % [Nx Ny]
    Npml = 1*[0, 15];  % [Nx_pml Ny_pml]
    mu0 = 4*pi*10^-7; mu_0 = mu0; mu = mu0;
    eps0 = 8.85*10^-12; eps_0 = eps0;
    [xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

    %% Set up the permittivity.
    eps_r = ones(N);
    eps_r(:,140:145) = -10+1i; %we need a nonzero imaginary part...to make the plasmon wavelength large?
    eps_r(240:250, 140:142) = 1;
    x = 1:N(1); y = 1:N(2);
    location_mask = zeros(N);
    % location_mask = speye(N);
    % location_mask = circshift(location_mask, -11);
    %location_mask(1:end, 40:end) = 0;

    location_mask(:,25) = 1;
    probe_loc = 25;
    figure; spy(location_mask);

    k = [1,0];
    Jz_plane = PlaneWaveSource(k, dL, location_mask,N);

    %% Set up the magnetic current source density.
    Mz = zeros(N);
    ind_src = [50,10];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
    Mz(ind_src(1), ind_src(2)) = 1;

    %% Run the solveTE Function
    [Hz, Ex, Ey, A, omega,b] = solveTE(L0, wvlen, xrange, yrange, eps_r, Jz_plane, Npml);
    visreal(Hz, xrange, yrange)
    figure
    visabs(Ex, xrange, yrange)
    figure
    visabs(Ey, xrange, yrange)
    figure;
%     moviereal(Hz, xrange, yrange)

    transmitted_probe = Hz(250,25);
    initial_probe = Hz(250,125);
    
    R = [R;(abs(transmitted_probe)/abs(initial_probe)).^2]
end
figure;
plot(R);