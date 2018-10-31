
%% test script for TM Ex Ey eigensolve
%% for anisotropy... the operator can start getting pretty dense...
% if that's the case, then a schur complement almost definitely makes sense
close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;
xrange = [-0.5,0.5];  % x boundaries in L0
yrange = [-1,1];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [200 200];  % [Nx Ny]
Npml = 1*[0 30];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);

%% Set up the permittivity.
%parfor wvlen = wvlen_scan
wvlen = 0.2;
delta = 5; delta2 = 50;
exx = ones(N);  exx(:, cy-delta:cy+delta) = 6;
eyy = ones(N);  eyy(:, cy-delta:cy+delta) = 2;
exy = zeros(N); exy(:, cy-delta:cy+delta) =0.2-0.01i;
eyx = zeros(N); eyx(:, cy-delta:cy+delta) = 0.2+0.01i;

% exx = ones(N);  exx(cx-delta2:cy+delta2, cy-delta:cy+delta) = 6;
% eyy = ones(N);  eyy(cx-delta2:cy+delta2, cy-delta:cy+delta) = 6;
% exy = zeros(N); exy(cx-delta2:cy+delta2, cy-delta:cy+delta) = 0.;
% eyx = zeros(N); eyx(cx-delta2:cy+delta2, cy-delta:cy+delta) = 0;

figure(); 
subplot(221);
imagesc(exx);
subplot(222);
imagesc(eyy);
subplot(223);
imagesc(abs(exy));
subplot(224);
imagesc(abs(eyx));
eps_tensor = {exx, exy; eyx, eyy};

%% eigensolve
neigs = 40;
[Hz_modes, Ex_modes, Ey_modes,~,A] = eigensolve_anisotropic_TM(L0, wvlen, xrange, ...
    yrange, eps_tensor, Npml, neigs);

for i = 1:neigs
   figure();
   subplot(121)
   visreal(Hz_modes{i}, xrange, yrange);
   subplot(122)
   visreal(Ey_modes{i}, xrange, yrange);
end

%% bloch eigensolve
% %% eigensolve
% Kx = pi/(2*L(1)); K_vec = [Kx, 0];
% neigs = 2;
% [Hz_modes, Ex_modes, Ey_modes,~] = ...
%     eigensolve_anisotropic_TM_bloch(L0, wvlen, xrange, ...
%     yrange, eps_tensor, Npml, neigs,K_vec);
% 
% for i = 1:neigs
%    figure();
%    subplot(121)
%    visreal(Ey_modes{i}, xrange, yrange);
%    subplot(122)
%    visreal(Ex_modes{i}, xrange, yrange);
% end
% 
% %% eigensolve with opposite moving K
% Kx = pi/(2*L(1)); K_vec = [-Kx, 0];
% neigs = 20;
% [Hz_modes, Ex_modes, Ey_modes, eigenvals] = ...
%     eigensolve_anisotropic_TM_bloch(L0, wvlen, xrange, ...
%     yrange, eps_tensor, Npml, neigs,K_vec);
% 
% for i = 1:neigs
%    figure();
%    subplot(121)
%    visreal(Ey_modes{i}, xrange, yrange);
%    subplot(122)
%    visreal(Ex_modes{i}, xrange, yrange);
% end
