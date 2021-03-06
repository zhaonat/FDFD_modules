%% test script for TM Ex Ey eigensolve

close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
mu0 = 4*pi*1e-7*L0;
eps0 = 8.854e-12*L0;
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
wvlen_guess = 0.5;

epsilon = ones(N);
x = 1:N(1);
y = 1:N(2);
[xx, yy] = meshgrid(x,y);
%epsilon((xx-cx).^2+(yy-cy).^2 < 50^2)=12;
epsilon(cx-40:cx+40, cy-40:cy+40) = -6-0.01i;
figure();
visabs(epsilon, xrange, yrange);
drawnow();

%% eigensolve
neigs = 40;
Kx = 0;
K_vec = [Kx,0];

[Hz, Ex, Ey, eigenvals,A] = eigensolve_TM(L0, wvlen_guess, xrange, ...
    yrange, epsilon, Npml, neigs, K_vec);

[Hzc, Exc, Eyc, eigenvalsc,Ac] = solveTE_BlochX(L0,...
    wvlen_guess, xrange, yrange, epsilon, Kx, Npml, neigs);

for i = 1:neigs
   figure();
   subplot(121)
   visreal(Hz{i}, xrange, yrange);
   subplot(122)
   visreal(Hzc{i}, xrange, yrange);
   title(eigenvals(i)/(2*pi*c0)*1e-6);

end

figure();
plot(real(eigenvals));
hold on;
plot(real(eigenvalsc));
figure();
plot(sort(abs(eigenvals)));
hold on;
plot(sort(abs(eigenvalsc)));

figure();
plot(abs(diag(A)))
hold on;
plot(abs(diag(Ac)))
