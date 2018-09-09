
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
Npml = 1*[40 40];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);

%% Set up the permittivity.
%parfor wvlen = wvlen_scan
x = 1:N(1);
y = 1:N(2);
[xx, yy] = meshgrid(x,y);
epsilon((xx-cx).^2+(yy-cy).^2 < 50^2)=12;
wvlen = 1; delta = 10;
exx = ones(N);  exx((xx-cx).^2+(yy-cy).^2 < 50^2) = 12;
eyy = ones(N);  eyy((xx-cx).^2+(yy-cy).^2 < 50^2) = 12;
exy = zeros(N); exy((xx-cx).^2+(yy-cy).^2 < 50^2) =6;
eyx = zeros(N); eyx((xx-cx).^2+(yy-cy).^2 < 50^2) = 2; %-1i;

% figure(); 
% subplot(221);
% imagesc(exx);
% subplot(222);
% imagesc(eyy);
% subplot(223);
% imagesc(exy);
% subplot(224);
% imagesc(eyx);
eps_tensor = {exx, exy; eyx, eyy};

%% eigensolve
Kx = pi/(2*L(1)); K_vec = [Kx, 0];
neigs = 2;
[Hz_modes, Ex_modes, Ey_modes, eigenvals] = ...
    eigensolve_anisotropic_TM_bloch(L0, wvlen, xrange, ...
    yrange, eps_tensor, Npml, neigs,K_vec);

for i = 1:neigs
   figure();
   subplot(121)
   moviereal(Hz_modes{i}, xrange, yrange);
   subplot(122)
   moviereal(Ex_modes{i}, xrange, yrange);
end

