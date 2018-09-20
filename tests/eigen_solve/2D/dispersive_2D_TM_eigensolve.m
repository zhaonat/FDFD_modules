%% test script for TM Ex Ey eigensolve

close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
mu0 = 4*pi*1e-7*L0;
eps0 = 8.854e-12*L0;
c0 = 3e8;
xrange = 0.2*[-1,1];  % x boundaries in L0
yrange = [-4,4];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [150 250];  % [Nx Ny]
Npml = 1*[0 20];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);

%% Set up the permittivity.
wvlen = 2;
omega = 2*pi*c0/(wvlen*L0);

epsilon = ones(N);
x = 1:N(1);
y = 1:N(2);
[xx, yy] = meshgrid(x,y);
%epsilon((xx-cx).^2+(yy-cy).^2 < 50^2)=12;
half_ny = 20;
epsilon(:,cy-half_ny:cy+half_ny) = 6;
epsilon(cx-10:cx+10, cy-half_ny:cy+half_ny) = -6-0.01i;
figure();
visabs(epsilon, xrange, yrange);
drawnow();

%% eigensolve
neigs = 30;

[Hz_modes, Ex_modes, Ey_modes, eigenvals] = ...
    eigensolve_TM_dispersive_Kx(L0, omega, xrange, yrange, epsilon, Npml, neigs);

for i = 1:neigs
    figure();
    Kx = eigenvals(i);
    visreal(Hz_modes{i}, xrange, yrange);
    title(strcat(num2str(i), ', ', num2str(real(Kx)/(2*pi)*diff(xrange))));
    
end

figure();
moviereal(Hz_modes{2}, xrange, yrange);