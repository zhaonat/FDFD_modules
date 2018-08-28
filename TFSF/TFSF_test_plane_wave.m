close all
clear
%% potential applications
% cavity mirrors: multi-frequency, lasers
% oled mirrors

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;
xrange = [-2 2];  % x boundaries in L0
yrange = [-2 2];  % y boundaries in L0
L = [diff(xrange), diff(yrange)]
N = [200 200];  % [Nx Ny]
Npml = 1*[0 40];  % [Nx_pml Ny_pml]
theta = pi/4;
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);

% for the source field
x = linspace(xrange(1), xrange(2), Nx);
y = linspace(yrange(1), yrange(2), Ny);
[X,Y] = meshgrid(x,y);
X = X.';
Y = Y.'

%% scattered field, total field mask: 0 = total, 1 = scattered;
Q = zeros(N); 
Q(:,end-80:end) = 1;      %scattered
figure();
imagesc(Q);
drawnow();
Q = diag(sparse(Q(:)));

%% Set up the permittivity.
spectra_r = [];
spectra_t = []; spectra_inc = [];
wvlen_scan = linspace(1.2,2, 20);
%parfor wvlen = wvlen_scan

wvlen = wvlen_scan(1);

k0 = 2*pi/(wvlen);
kx = k0*cos(theta);
ky = k0*sin(theta);;
thickness = 1;

epsilon= ones(Nx,Ny);

%     figure();
%     imagesc(abs(eps));
%     drawnow()
%% Set up the magnetic current source density.
Mz = zeros(N);

%% source location
ind_src = [60,220];  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(:, ind_src(2)) = 1/Nx;



%% Solve TE (photonic), TM (practical), depends on convention equations.
tic
[Ez, Hx, Hy, A, b] = solveTM(L0, wvlen, xrange, yrange, epsilon, Mz, Npml);
toc

%% source field calculation
fsrc = A\b; %this has no particular localization...

Fsrc = (sparse(fsrc(:)));

bsrc = (Q*A-A*Q)*Fsrc; %Fsrc is not the same as b... which is interesting...
tic
Ez_tfsf = A\bsrc;
Ez_tfsf = reshape(Ez_tfsf,Nx,Ny);
toc
figure()
subplot(121)
moviereal(Ez_tfsf, xrange, yrange)


%legend('ref','tran')



