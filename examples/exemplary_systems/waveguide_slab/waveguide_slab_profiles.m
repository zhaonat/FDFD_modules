%% Set up the domain parameters.
close all

L0 = 1e-9;  % length unit: microns
eps_0 = 8.85*10^-12*L0; 
mu_0 = 4*pi*10^-7*L0; 
c0 = 1/(eps_0*mu_0)^.5;
xrange = [-400 400];  % x boundaries in L0
yrange = [-2000, 2000];  % y boundaries in L0
N = [160 320];  % [Nx Ny]
Npml = 1*[0 15];  % [Nx_pml Ny_pml]
Nx = N(1); Ny = N(2);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
eps_r = ones(N);
epsilon = 12; divet_width = 100; divet_height = 500; slab_width = 1000;

eps_divet = slab_extrusion(eps_r,dL, epsilon, ...
    divet_width, divet_height, slab_width);
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
    yrange, eps_r, Npml, num_eigs, K_vec);

% mode visualization
for i  = 1:num_eigs
    mode = Ez_modes{i};
    ratio = energy_distribution(mode, structure_mask);
    if(ratio > 0.2)
        figure;
        subplot(121)
        visreal(mode, xrange, yrange);
        subplot(122);
        plot(real(mode(floor(N(1)/2)-2: floor(N(1)/2)+2,:)).');
        title(ratio)
    end
end
