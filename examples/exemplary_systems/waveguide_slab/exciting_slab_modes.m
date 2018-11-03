%% Set up the domain parameters.
close all
clear

L0 = 1e-9;  % length unit: microns
eps_0 = 8.85*10^-12*L0; 
mu_0 = 4*pi*10^-7*L0; 
c0 = 1/(eps_0*mu_0)^.5;
xrange = [-400 400];  % x boundaries in L0
yrange = [-2000, 2000];  % y boundaries in L0
N = [160 320];  % [Nx Ny]
Npml = 1*[0 15];  % [Nx_pml Ny_pml]
Nx = N(1); Ny = N(2);
[xrange, yrange1, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
eps_r = ones(N);
epsilon = 12; slab_width = 1000;
eps_slab = dielectricSlab(eps_r, dL, slab_width, epsilon);

eps_r = eps_slab;

wvlen = 2000;
omega_guess = 2*pi*c0/wvlen;
Kx = 0.5*(pi/diff(xrange));
figure; 
imagesc(eps_r)
K_vec = [Kx, 0];
num_eigs = 50;

structure_mask = (eps_r == epsilon);
[Ez_modes, Hx_modes, Hy_modes, eigenvals] = eigensolve_TE(L0, wvlen, xrange, ...
    yrange1, eps_r, Npml, num_eigs, K_vec);

for i  = 1:num_eigs
    mode_ex = Ez_modes{i};
    ratio = energy_distribution(mode_ex, structure_mask);
    if(ratio > 0.2)
        figure;
        subplot(121)
        visreal(mode_ex, xrange, yrange1);
        subplot(122);
        plot(real(mode_ex(floor(N(1)/2)-2: floor(N(1)/2)+2,:)).');
        title(i)
    end
end

%% Try to excite it by placing a mode as a source;
mode_index = 43;
figure(); visreal(Ez_modes{mode_index }, xrange, yrange1);
figure(); plot(real(Ez_modes{mode_index}(floor(N(1)/2),:)));
Npml_sim = 1*[0, 15];  % [Nx_pml Ny_pml]

xrange2 = [-2000 2000];  % x boundaries in L0
yrange = [-2000, 2000];  % y boundaries in L0
N2 = [600 320];  % [Nx Ny]

[xrange2, yrange2, N2, dL2, Lpml] = domain_with_pml(xrange2, yrange, N2, Npml);  % domain is expanded to include PML
eps_r = ones(N2);
epsilon = 12; 
eps_slab = dielectricSlab(eps_r, dL2, slab_width, epsilon);

eps_r = eps_slab;
Mz = zeros(N2);

Mz(80,:) = Ez_modes{mode_index}(80,:);
figure(); visreal(Mz, xrange2, yrange2);

wvlen = 1985;
[Hz, Ex, Ey, A] = ...
    solveTM(L0, wvlen, xrange2, yrange2, eps_r, Mz, Npml_sim);
figure();moviereal(Hz, xrange2, yrange2);
