%% 1D FDFD simulation
clear
close all

%% === ESSENTIAL SIMULATION PARAMETERS=================
L0 = 1e-6;  % length unit: microns
wvlen = 1.55;  % wavelength in L0
iterSchur = []; iterUnRed = [];
maxCellNumber = 1;
SingleCellSize = 50;
epsilon = 1;
Sx = SingleCellSize;
Nx = SingleCellSize;

%% NPML construction
Npml = [0,0];

%% epsilon dielectric construction;
eps_r = ones(1, Nx);
eps_r(15:35) = 1;

%% Simulation Domain
zrange = [-1 1];  % x boundaries in L0
wvlen = 1;    
L0 = 10^-6;
dz = (zrange(2)-zrange(1))/Nx;
dL = [dz];

%% establish Source
Mz = zeros(Nx,1);
ind_src = [ceil(Nx/2) ceil(Nx/2)];%ceil(AN/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), 1) = 1;

%% 1D simulation
[A, omega,b, Dxf, Dxb, Tep] = ...
    solveFDFD_1D(wvlen, zrange, dL, eps_r, Mz, Npml, L0);

%% diagnose A
figure; spy(A);

%% Test Solve
figure;
plot(abs(A\b));

%% Schur reordering
[SymA, SymB, Q] = Reordering1D(A,b);

%% Perform Schur complement
hpart = 2; vpart = 2;
[Aschur, bmod, App, Apv, Avp, Avv, bp, bv] = ...
    FourblockSchur(SymA,SymB, hpart, vpart)

xschur = Aschur\bmod;
xtol = SymA\SymB;

