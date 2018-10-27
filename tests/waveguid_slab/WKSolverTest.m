%% Set up the domain parameters.
close all

L0 = 1e-9;  % length unit: microns
eps_0 = 8.85*10^-12*L0; 
mu_0 = 4*pi*10^-7*L0; 
c0 = 1/(eps_0*mu_0)^.5;
xrange = [-400 400];  % x boundaries in L0
yrange = [-2000, 2000];  % y boundaries in L0
N = [40 160];  % [Nx Ny]
Npml = 1*[0 15];  % [Nx_pml Ny_pml]
Nx = N(1); Ny = N(2);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
eps_r = ones(N);
epsilon = 12; divet_width = 100; divet_height = 500; slab_width = 500;

eps_divet = slab_extrusion(eps_r,dL, epsilon, ...
    divet_width, divet_height, slab_width);
eps_slab = dielectricSlab(eps_r, dL, slab_width, epsilon);

eps_r = eps_divet;

wvlen = 2000;
omega_guess = 2*pi*c0/wvlen;
K = 0.5*(pi*L0);
figure; 
imagesc(eps_r)

structure_mask = (eps_r == epsilon);
A = TEWKSolver(L0, xrange, yrange, eps_r, Npml,K, omega_guess);

size(A)
num_eigs = 50;
[U,S,V] = eigs(A, num_eigs, 'sm');

% mode visualization
for i  = 1:num_eigs
mode = reshape(U(:,i),N(1), N(2));

    ratio = energy_distribution(mode, structure_mask);
    if(ratio > 0.2)
        figure;
        visreal(mode, xrange, yrange);
        title(ratio)
    end
end
