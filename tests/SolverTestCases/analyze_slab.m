% EE256:Numerical Electromagnetics, Spring 2013-2014
% Homework 4
% Problem 2: Revisiting the transmission through a slab
% Script by Chunliang Zheng

clear all; close all; clc;

%% Set up the domain parameters.
L0 = 1e-9;  % length unit: nm
xrange = [-300 1000];  % x boundaries in L0
yrange = [0 500];  % y boundaries in L0
N = [260 1];  % [Nx Ny]
Npml = [10 0];  % [Nx_pml, Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
d = 700;  % slab thickness
eps_slab = 12;  % slab permittivity
within_slab = @(x,y) x > 0 & x < d;  % function handle that returns true if x is in slab
eps_with_slab = ones(N);
eps_with_slab = assign_val(eps_with_slab, xrange, yrange, within_slab, eps_slab);

eps_no_slab = ones(N);

%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = round((Lpml(1) + 100) / dL(1));
Mz(ind_src, :) = 1;


%% Print out SWR and the upper bound of rho.
wvlen=800;
[Hz, Ex, Ey, A, omega] = solveTE(L0, wvlen, xrange, yrange, eps_no_slab, Mz, Npml);
figure;
visabs(Hz, xrange, yrange);
figure;
moviereal(Hz, xrange, yrange)

swr = max(abs(Hz(Npml(1)+1:N-Npml(1))))/min(abs(Hz(Npml(1)+1:N-Npml(1)))); 
rho_max = (swr-1)/2;
disp(['SWR = ', num2str(swr)]);
disp(['rho <= ', num2str(rho_max)])


%% Plot the measured and theoretical transmittance spectra.
n_slab=sqrt(eps_slab);
i=1;
for wvlen=730:2:870
    [Hz_no_slab, Ex_no_slab, Ey_no_slab, A_no_slab, omega_no_slab] = solveTE(L0, wvlen, xrange, yrange, eps_no_slab, Mz, Npml);
    [Hz_with_slab, Ex_with_slab, Ey_with_slab, A, omega_with_slab] = solveTE(L0, wvlen, xrange, yrange, eps_with_slab, Mz, Npml);
    Ts_measure(i)=abs(Hz_with_slab(251)/Hz_no_slab(251))^2;
    Ts_theory(i)=abs(4*n_slab/((n_slab-1)^2-exp(1i*4*pi*n_slab*d/wvlen)*(n_slab+1)^2))^2;
    i=i+1;
end

wvlens=730:2:870;
figure;
plot(wvlens, Ts_theory, 'b', wvlens, Ts_measure, 'ro')
legend('Theory', 'FDFD', 'Location', 'NW');


%% Grading breakdown: (total 2 points)

% You get 1 point if you calculate swr and rho_max correctly
% You get 1 point if you plot the measured transmission and theoretical transmission correctly
  
