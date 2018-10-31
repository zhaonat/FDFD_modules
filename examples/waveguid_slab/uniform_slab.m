%% TYPICAL bAND sTRUCTURE OF A GUIDED RESONANCE SLAB

close all
clear
%% Set up the domain parameters.
L0 = 1e-9;  % length unit: nanometers
eps_0 = 8.85*10^-12*L0;
mu_0 = 4*pi*10^-7*L0; 
c0 = 1/(eps_0*mu_0)^.5;
xrange = [-200 200];  % x boundaries in L0
yrange = [-2000, 2000];  % y boundaries in L0
N = [40 160];  % [Nx Ny]
Npml = 1*[0 15];  % [Nx_pml Ny_pml]
epsilon = 40;
Nx = N(1); Ny = N(2);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
eps_r = ones(N);
slab_width = 600;
eps_slab = dielectricSlab(eps_r, dL, slab_width, epsilon);
divetSize = 8; width = 20;
%eps_divet = WaveGuideDivet(eps_r,epsilon, divetSize, width)
inner_rad = 0; outer_rad = 40
%eps_divet = RingResonator(eps_r, inner_rad, outer_rad, epsilon);
eps_r = eps_slab;
wvlen = 500;
omega_guess = 2*pi*c0/wvlen;
K = 0.5*(pi*L0);;
figure; 
imagesc(eps_r)

structure_mask = (eps_r == epsilon);
figure; 
for K = linspace(0, pi/(diff(xrange)), 200)
    frequencies = [];
    A = TEWKSolver(L0, xrange, yrange, eps_r, Npml,K, omega_guess);
    num_eigs = 10;
    [U,S,V] = eigs(A, num_eigs, 'sm');
    s = diag(S);
    % mode visualization
    for i  = 1:10
        mode = reshape(U(:,i),N(1), N(2));
       
        ratio = energy_distribution(mode, structure_mask)
        if(ratio> 0.15)
            frequencies = [frequencies, s(i)];
        end
    end
    scatter(repmat(K, length(frequencies), 1), abs(frequencies).^.5, '.b');
    hold on;

end

%% plot the light line
hold on;
Kspace = linspace(0, pi/(diff(xrange)), 100);
plot(Kspace, c0*Kspace);
%% plot the medium light line
hold on
plot(Kspace, (c0/sqrt(epsilon))*Kspace)