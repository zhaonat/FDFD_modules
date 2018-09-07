
%% test script for TM Ex Ey eigensolve
close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;
xrange = [-0.2,0.2];  % x boundaries in L0
yrange = [-1,1];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [500 500];  % [Nx Ny]
Npml = 1*[0 40];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);

%% Set up the permittivity.
%parfor wvlen = wvlen_scan
wvlen = 1;
exx = ones(N);  exx(:, cy-50:cy+50) = 12;
eyy = ones(N);  eyy(:, cy-50:cy+50) = 12;
exy = zeros(N); exy(:, cy-50:cy+50) = 0;
eyx = zeros(N); eyx(:, cy-50:cy+50) = 0;

figure(); 
subplot(221);
imagesc(exx);
subplot(222);
imagesc(eyy);
subplot(223);
imagesc(exy);
subplot(224);
imagesc(eyx);
eps_tensor = {exx, exy; eyx, eyy};

%% eigensolve

neigs = 20;
[Hz_modes, Ex_modes, Ey_modes, eigenvals] = eigensolve_anisotropic_TM(L0, wvlen, xrange, ...
    yrange, eps_tensor, Npml, neigs);

for i = 1:neigs
   figure();
   subplot(121)
   visreal(Hz_modes{i}, xrange, yrange);
   subplot(122)
   visreal(Ex_modes{i}, xrange, yrange);
end
