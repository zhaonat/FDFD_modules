%% test script for TM Ex Ey eigensolve

close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;
xrange = [-1,1];  % x boundaries in L0
yrange = [-1,1];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [200 200];  % [Nx Ny]
Npml = 1*[20 20];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);

%% Set up the permittivity.
wvlen = 1;

epsilon = ones(N);
x = 1:N(1);
y = 1:N(2);
[xx, yy] = meshgrid(x,y);
%epsilon((xx-cx).^2+(yy-cy).^2 < 50^2)=12;
epsilon(cx-40:cx+40, cy-40:cy+40) = -3-0.01i;
%epsilon(cx-40:cx+40, cy+45:cy+105) = -3-0.01i;

%epsilon(cx+55:cx+105, cy+55:cy+105) = -3-0.01i;
figure();
visabs(epsilon, xrange, yrange);
drawnow();

%% source
Mz = zeros(N);
Mz(75,cy) = 1;


[Hz, Ex, Ey] = solveTE(L0, wvlen, xrange, yrange, epsilon, Mz, Npml);


% Kx = 0;
% [Hzc, Exc, Eyc, eigenvalsc] = solveTE_BlochX(L0,...
%     wvlen_guess, xrange, yrange, epsilon, Kx, Npml, neigs);

figure()
visreal(Hz, xrange, yrange);
figure(); 
moviereal(Hz, xrange, yrange);
