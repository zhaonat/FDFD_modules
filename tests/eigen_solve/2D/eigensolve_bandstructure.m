%% test script for TM Ex Ey eigensolve

close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
mu0= 4*pi*1e-7;
c0 = 3e8;
xrange = 0.1*[-1,1];  % x boundaries in L0
yrange = [-1,1];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [200 200];  % [Nx Ny]
Npml = 1*[0 20];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);

a = L(1)+Lpml(1);
%% Set up the permittivity.
wvlen_guess = 1;
epsilon = ones(N);
x = 1:N(1);
y = 1:N(2);
[xx, yy] = meshgrid(x,y);
epsilon(:, cy-20:cy+20)=12;
figure();
imagesc(epsilon);
drawnow();

%% eigensolve
neigs = 4;
omega_eigs = [];
omega_eigs_check = [];
kx_scan = linspace(0, pi/a, 2);
for Kx = kx_scan
    Kx
    K_vec = [Kx,0];

    [Hz, Ex, Ey, eigenvals,A] = eigensolve_TM(L0, wvlen_guess, xrange, ...
        yrange, epsilon, Npml, neigs, K_vec);

    [Hzc, Exc, Eyc, eigenvalsc,Ac] = solveTE_BlochX(L0,...
        wvlen_guess, xrange, yrange, epsilon, Kx, Npml, neigs);
    omega_eigs = [omega_eigs, eigenvals];
    omega_eigs_check = [omega_eigs_check, eigenvalsc];
end

figure();
plot(kx_scan, abs(omega_eigs));


figure();
plot(kx_scan, abs(omega_eigs_check));

