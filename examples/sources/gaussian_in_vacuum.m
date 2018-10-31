close all
clear
%% potential applications
% cavity mirrors: multi-frequency, lasers

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;
wvlen = 1;
xrange = [-3 3];  % x boundaries in L0
yrange = [-2 2];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [300 200];  % [Nx Ny]
Npml = 1*[40 40];  % [Nx_pml Ny_pml]
y1 = 1.5;
y2 = 2.5;

[xrange, yrange, N, dL, Lpml] =domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);

%% Set up the permittivity
eps = ones(N);

figure(); visabs(eps, xrange, yrange);
drawnow();
%% Set up the magnetic current source density.
sigma = 0.2; n_source = 150;
Kx = 2*pi/wvlen;
Mz = gaussian_grid(N,Npml,xrange, yrange, sigma, Kx);

Mz = zeros(N);
Nx = 240;
xspace = linspace(-4,4,Nx);
sigma = 1;
line_source = gaussian_line(xspace, Nx, sigma)
Mz(350,140-120:140+120-1) = line_source;
figure();
visreal(Mz,xrange,yrange);
drawnow();

%% set up electric current source density
Jx = zeros(N);
Jy = zeros(N);

%% Solve TE (photonic), TM (practical), depends on convention equations.
tic
[Hz, Ex, Ey] = solveTE(L0, wvlen, xrange, yrange, eps, Mz, Npml);
toc

figure();
visabs(Hz, xrange, yrange)

figure();
visreal(Hz, xrange, yrange);

figure();
moviereal(Hz,xrange,yrange);



