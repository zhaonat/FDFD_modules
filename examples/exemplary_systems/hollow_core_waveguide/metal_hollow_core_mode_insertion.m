%% test script for TM Ex Ey eigensolve
%% RESULTS: 11/2/2018; shows that you can take the mode source and correctly
% excite it in a driven simulation

close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
mu0 = 4*pi*1e-7*L0;
eps0 = 8.854e-12*L0;
c0 = 3e8;
lattice_constant = 0.2;
num_cells = 15;

xrange = num_cells*(lattice_constant/2)*[-1,1];  % x boundaries in L0
yrange = 3*[-1,1];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [500 300];  % [Nx Ny]
Npml = 1*[0 15];  % [Nx_pml Ny_pml]

%[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
dL = L./N;
Lpml = [0,0];
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);
wvlen = 1.8;
fill_factor = 0.2;

%% Set up the permittivity.
omega = 2*pi*c0/(wvlen*L0);
omega_p =  0.72*pi*1e15;
gamma = 5.5e12;
epsilon_metal = 1-omega_p^2/(omega^2-1i*gamma*omega);
epsilon_array = [16, epsilon_metal];

% our convention on the dimensions is not good enough...
ycentercore = 0;
hollow_core_half_thickness = 1.5;
wall_thickness = 1;
y1 =  ycentercore-hollow_core_half_thickness-wall_thickness/2 ;
y2 = ycentercore+hollow_core_half_thickness+wall_thickness/2 ;

air_core = air_core_structure(xrange, yrange,N, Lpml);
wall_properties = {num_cells, lattice_constant, ...
                            wall_thickness, epsilon_array, fill_factor};
wall_properties2 = {num_cells, lattice_constant, ...
                            wall_thickness, epsilon_array, fill_factor};
                        
y_centers = [y1, y2]; %these are the centers of the WALL
air_core. build_multi_core(y_centers, {wall_properties, wall_properties2})

figure();
visreal(air_core.epsilon, xrange, yrange);
drawnow();

%% filtering bounds
xlim = xrange;
ylim = [ycentercore-hollow_core_half_thickness, ycentercore+hollow_core_half_thickness];
epsilon = air_core.epsilon;

%% eigensolve
neigs = 120;
kx_guess = 0*pi/L(1);
[Hz_modes, Ex_modes, Ey_modes, kx_eigs] = ...
   eigensolve_TM_dispersive_Kx(L0, omega, xrange, yrange, epsilon, Npml, neigs, kx_guess);
[filtered_modes, filtered_k] = ...
    mode_filtering(Hz_modes, kx_eigs, epsilon, xlim, ylim, L, Npml);

% K_vec = [0.5*pi/L(1),0];
% [Hz_modes, Ex_modes, Ey_modes, eigenvals, A] = eigensolve_TM(L0, wvlen, xrange, ...
%     yrange, epsilon, Npml, neigs, K_vec);
% 
% [filtered_modes, filtered_k] = ...
%     mode_filtering(Hz_modes, eigenvals, epsilon, xlim, ylim, L, Npml);

for i = 1:length(filtered_k)
    figure();
    Kx = filtered_k(i);%small Kx, but these modes should not be asymmetric
    visreal(filtered_modes{i}, air_core.xrange, air_core.yrange);
    title(strcat(num2str(i), ', ', num2str(imag(Kx)/(2*pi)*diff(xrange))));  
end
drawnow();


%% EXCITING USING MODE SOURCE
index =18;
figure(); 
subplot(121)
visabs(filtered_modes{index}, air_core.xrange, air_core.yrange);
subplot(122)
visreal(filtered_modes{index}, air_core.xrange, air_core.yrange)
source = filtered_modes{index};

%% run the final simulation
%Mz(floor(N(1)/2), floor(N(2)/2)) = 1;
lattice_constant = 0.2;
num_cells = 15;
xrange = num_cells*(lattice_constant/2)*[-1,1];  % x boundaries in L0
L = [diff(xrange), diff(yrange)];
N_sim = [500 300];  % [Nx Ny]
dL = L./N_sim;

long_air_core = air_core_structure(xrange, yrange,N_sim, Lpml);
wall_properties = {num_cells, lattice_constant, ...
                            wall_thickness, epsilon_array, fill_factor};
wall_properties2 = {num_cells, lattice_constant, ...
                            wall_thickness, epsilon_array, fill_factor};
                        
y_centers = [y1, y2]; %these are the centers of the WALL
long_air_core. build_multi_core(y_centers, {wall_properties, wall_properties2})

figure(); 
subplot(121)
visreal(long_air_core.epsilon, long_air_core.xrange, long_air_core.yrange);
subplot(122)
visreal(air_core.epsilon, air_core.xrange, air_core.yrange);

Mz = zeros(N_sim);

Npml2 = [0,15];
Mz(100,:) = (source(floor(N(1)/2),:));
wvlen = 2
[Hz, Ex, Ey] = ...
    solveTE(L0, wvlen, xrange, yrange, long_air_core.epsilon, Mz, Npml2);

figure();
visreal(Hz, xrange, yrange);
figure();
visabs(Hz, xrange, yrange);

figure(); 
subplot(121)
visreal(Mz, xrange, yrange);
line(xrange, [y1, y1]+wall_thickness/2)
line(xrange, [y2, y2]-wall_thickness/2)

subplot(122)
plot(real(source(70,:))/max(abs(source(70,:))));
hold on;
plot(real(Hz(70,:))/max(abs(Hz(70,:))));
