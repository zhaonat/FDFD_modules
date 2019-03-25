close all
clear
%% potential applications
% cavity mirrors: multi-frequency, lasers
% oled mirrors

%% NOTE THAT PMLs cannot handle glancing incidence!!! 
%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;
xrange = [-1 1];  % x boundaries in L0
yrange = [-2 2];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [202 300];  % [Nx Ny]
Npml = 1*[0,15];  % [Nx_pml Ny_pml]
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

epsilon = ones(N);
%% Set up the magnetic current source density.
wvlen = 1;

Mz = zeros(N);
ind_src = [1, N(2)-Npml(2)-40];  % (i,j) indices of the center cell; Nx, Ny should be odd
degrees = pi/180;
theta = 70*degrees;

kx = 2*pi/wvlen*sin(theta);
ky = 2*pi/wvlen*cos(theta);
Kvec = [kx, ky];

x_vec = linspace(xrange(1),xrange(2), N(1));
Mz(:, ind_src(2)) = exp(-1i*(kx*x_vec)); %this sets kx = 0; or theta = 0;

Mz(:, ind_src(2)+1) = Mz(:, ind_src(2)).*exp(-1i*(ky * dL(2)) - 1i*pi);
    
figure();
visreal(Mz, xrange, yrange);

%% Solve TE (photonic), TM (practical), depends on convention equations.
tic
[Hz, Ex, Ey] = ...
    solveTE_bloch(L0, wvlen, xrange, yrange, epsilon, Mz, Npml, Kvec);
toc


%% figure
figure();
moviereal(Hz, xrange, yrange);