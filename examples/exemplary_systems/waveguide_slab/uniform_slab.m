%% TYPICAL bAND sTRUCTURE OF A GUIDED RESONANCE SLAB

close all
clear
%% Set up the domain parameters.
L0 = 1e-9;  % length unit: nanometers
eps_0 = 8.85*10^-12*L0;
mu_0 = 4*pi*10^-7*L0; 
c0 = 1/(eps_0*mu_0)^.5;
xrange = [-100 100];  % x boundaries in L0
yrange = [-2000, 2000];  % y boundaries in L0
N = [40 160];  % [Nx Ny]
Npml = 1*[0 15];  % [Nx_pml Ny_pml]
epsilon = 2; 
Nx = N(1); Ny = N(2);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
eps_r = ones(N);
slab_width = 1000;
eps_slab = dielectricSlab(eps_r, dL, slab_width, epsilon);
divetSize = 8; width = 20;
%eps_divet = WaveGuideDivet(eps_r,epsilon, divetSize, width)
inner_rad = 0; outer_rad = 40
%eps_divet = RingResonator(eps_r, inner_rad, outer_rad, epsilon);
eps_r = eps_slab;
wvlen = 500;
omega_guess = 2*pi*c0/wvlen;
figure; 
visreal(eps_r, xrange, yrange);

structure_mask = (eps_r == epsilon);
figure; 
num_eigs = 30;
for K = linspace(-pi/(2*diff(xrange)), pi/(2*diff(xrange)), 80)
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
    scatter(repmat(K, length(frequencies), 1), frequencies, '.b');
    hold on;
    drawnow();

end

%% plot the light line
hold on;
Kspace = linspace(0, pi/(diff(xrange)), 100);
plot(Kspace, c0*Kspace);
%% plot the medium light line
hold on
plot(Kspace, (c0/sqrt(epsilon))*Kspace)