close all

%% specify a simple FDFD
c0 = 3*10^8;
L0 = 1e-6;  % length unit: microns
wvlen = 3.0*L0;  % wavelength in L0
xrange = [-5 5]*L0;  % x boundaries in L0
yrange = [-5 5]*L0;  % y boundaries in L0
Npml = 1*[15 15];  % [Nx_pml Ny_pml]
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% TFSF Partition
Nx = 100; Ny = 100;
N = [Nx,Ny];

Q = ones(Nx,Ny); %0 indexes the total field
Q(25:75, 25:75) = 0; %1 indexes the scattered field, so we isolate only scattered field
Q(45:55, 45:55) = 1;
% convert into diagonal matrix sparse(
M = prod(N);
Qr = diag(sparse(Q(:)));

%% Set up the permittivity.
eps_r = ones(N);

%% Set up the 1agnetic current source density.
Jz = zeros(N);
%must be entirely inside the TFSF interface
ind_src = [50,50];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Jz(ind_src(1), ind_src(2)) = 1;

%%solve for the source field in vacuum
Npml0 = [0,0];
[Ez, Hx, Hy, A0, b] = solveTM(L0, wvlen, xrange, yrange, eps_r, Jz, Npml);
fsrc = A0\b; %vacuum field distribution of the point source


Fsrc = reshape(fsrc, Nx,Ny);
figure; 
visreal(Fsrc, xrange, yrange);

%% ADD SCATTERER
eps_r(20:40, 20:40) = 12;
[Eza, Hxa, Hya, A,b] = solveTM(L0, wvlen, xrange, yrange, eps_r, Jz, Npml);

bprime = (Qr*A - A*Qr)*fsrc;

%% solve TFSF problem
x = A\bprime; %% the tfsf problem isolates ONLY  the scattered fields...?

e = reshape(x, Nx,Ny);
visreal(1i*e, xrange, yrange);

figure;
visreal(1i*Eza, xrange, yrange);


