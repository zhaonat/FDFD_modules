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
epsilon = 12; divet_width = 100; divet_height=100; slab_width = 500;

eps_divet = slab_extrusion(eps_r,dL, epsilon, ...
    divet_width, divet_height, slab_width);
eps_r = eps_divet;
wvlen = 200;
omega_guess = 2*pi*c0/wvlen;
K = 0.5*(pi*L0); %% K which is roughly pi/a, has to be scaled by L0
figure; 
imagesc(eps_r)

structure_mask = (eps_r == epsilon);
figure; 
BS = [];
for K = linspace(0, pi/(diff(xrange)), 100)
    frequencies = [];
    A = TEWKSolver(L0, xrange, yrange, eps_r, Npml,K, omega_guess);
    num_eigs = 10;
    [U,S,V] = eigs(A, num_eigs, 'sm');
    s = diag(S);
    % mode visualization
    for i  = 1:10
        mode = reshape(U(:,i),N(1), N(2));
        ratio = energy_distribution(mode, structure_mask)
        if(ratio> 0.1)
            frequencies = [frequencies, s(i)];
        end
    end
    scatter(repmat(K, length(frequencies), 1), abs(frequencies).^(.5), '.b');
    drawnow();
    hold on;
    BS = [BS, abs(frequencies).^(.5)];
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
