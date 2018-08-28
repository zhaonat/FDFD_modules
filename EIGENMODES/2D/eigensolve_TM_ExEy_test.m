

%% test script for TM Ex Ey eigensolve
close all
clear
%% potential applications
% cavity mirrors: multi-frequency, lasers
% oled mirrors
% it appears that Jx and Jy dipoles 

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;
xrange = [-0.2,0.2];  % x boundaries in L0
yrange = [-1,1];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [40 40];  % [Nx Ny]
Npml = 1*[0 0];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);
%% Set up the permittivity.

%parfor wvlen = wvlen_scan
reflected_power = [];
wvlen = 1.7;

exx = ones(N); exx(:, cy-10:cy+10) = 12;
eyy = ones(N); eyy(:, cy-10:cy+10) = 12;
exy = ones(N);exy(:, cy-10:cy+10) = 2;
eyx = ones(N);eyx(:, cy-10:cy+10) = -2;

eps_tensor = {exx, exy; eyx, eyy};
%% eigensolve
neigs = 10;

% [eigenvals, eigenmodes,A] = eigensolve_TM_Hz(L0, wvlen, xrange, ...
%     yrange, eps_tensor, Npml, neigs);
Kx = pi/L(1)/2;
[eigenvals, eigenmodes,A] = eigensolve_anisotropic_TM_Hz_bloch(L0, wvlen, xrange, ...
    yrange, eps_tensor, Npml, Kx, neigs, L);
% 
% figure(); visreal(eigenmodes{3,1}, xrange, yrange);
% 
% for i = 1:neigs
%    figure();
%    visreal(eigenmodes{i,2}, xrange, yrange);
% end

Kx = -pi/L(1)/2;
[eigenvals2, eigenmodes2,A2] = eigensolve_anisotropic_TM_Hz_bloch(L0, wvlen, xrange, ...
    yrange, eps_tensor, Npml, Kx, neigs, L);
for i = 1:neigs
   figure();
   subplot(121)
   moviereal(eigenmodes2{i,2}, xrange, yrange);
   subplot(122)
   moviereal(eigenmodes{i,2}, xrange, yrange);

end
