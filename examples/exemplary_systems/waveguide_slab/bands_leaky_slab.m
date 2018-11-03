%%GUIDED RESONANCE SLAB
close all
clear
%% Set up the domain parameters.
L0 = 1e-9;  % length unit: microns
eps_0 = 8.85*10^-12*L0; eps0 = 8.85*106-12;
mu_0 = 4*pi*10^-7*L0; mu0 = 4*pi*10^-7;
c0 = 1/(eps_0*mu_0)^.5;
xrange = [-150 150];  % x boundaries in L0
yrange = [-2000, 2000];  % y boundaries in L0
N = [80 160];  % [Nx Ny]
Npml = 1*[0 15];  % [Nx_pml Ny_pml]
Nx = N(1); Ny = N(2);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
eps_r = ones(N);
epsilon = 2; divet_width = 100; divet_height=100; slab_width = 1000;

eps_divet = slab_extrusion(eps_r,dL, epsilon, ...
    divet_width, divet_height, slab_width);
eps_r = eps_divet;
wvlen = 2000;
omega_guess = 2*pi*c0/wvlen;
K = 0.5*(pi*L0); %% K which is roughly pi/a, has to be scaled by L0
figure; 
imagesc(eps_r)

structure_mask = (eps_r == epsilon);
figure; 
num_eigs = 40;
BS = [];
for K = linspace(0, pi/(2*diff(xrange)), 50)
    K_vec = [K,0];
    frequencies = [];
    [Ez_modes, Hx_modes, Hy_modes, eigenvals] = eigensolve_TE(L0, wvlen, xrange, ...
    yrange, eps_r, Npml, num_eigs, K_vec);
    % mode visualization
    for i  = 1:num_eigs
        mode = Ez_modes{i};
       
        ratio = energy_distribution(mode, structure_mask); % energy_inside/energy_total
        if(ratio> 0.1)
            frequencies = [frequencies, eigenvals(i)];
        end
    end
    scatter(repmat(K, length(frequencies), 1), real(frequencies), '.b');
    drawnow();
    hold on;
    BS = [BS, real(frequencies)];
    %hold on;

end
%% plot the light line
hold on;
Kspace = linspace(0, pi/(diff(xrange)), 100);
plot(Kspace, c0*Kspace);
%% plot the medium light line
hold on
a = diff(xrange);
plot(Kspace, (c0/sqrt(epsilon))*Kspace)
plot(Kspace, (c0/sqrt(epsilon))*(Kspace+2*pi/a))
plot(Kspace, (c0/sqrt(epsilon))*(Kspace+4*pi/a))
