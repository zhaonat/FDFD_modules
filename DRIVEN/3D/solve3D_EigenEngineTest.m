%% Convergence History Test for the Wonsoek Accelerator solve3d_EigenEngine
c0 = 3*10^8;
L0 = 1e-6;  % length unit: microns
wvlen = 35.0*L0;  % wavelength in L0
xrange = [-5 5]*L0;  % x boundaries in L0
yrange = [-5 5]*L0;  % y boundaries in L0
zrange = [-5 5]*L0;
i = 10;
Nx = i; Ny = i; Nz = i;
N = [Nx Ny Nz];  % [Nx Ny]
M = N(1)*N(2)*N(3);
eps_r = ones(N);

%% Set up the current source density.
Mz = zeros(N); My = Mz; Mx = Mz;
ind_src = [2 2 2];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2), ind_src(3)) = 1;
JCurrentVector = [Mx; My; Mz];

%% PML boundary condition not yet implemented in the 3D case
Npml = [0 0 0];  % [Nx_pml Ny_pml]
mu0 = 4*pi*10^-7; mu_0 = mu0; mu = mu0;
eps0 = 8.85*10^-12; eps_0 = eps0;
tol = 1e-6;
figure;
for s = -2:1:2
   s
   residuals = [];
   [Ex, Ey, Ez, A,b] = solve3D_EigenEngine(L0, wvlen, xrange, yrange, zrange, eps_r, JCurrentVector, Npml,s);
   
   for maxit = 1:100
       [x1, flag, iter1, Miter, QlPiter, relres1, relAres,Anorm, Acond1, xnorm, axnorm, resvec1] = minresQLP(A,b,tol, maxit);    
       residuals = [residuals relres1];
   end
   
   figure;
   plot(log(residuals));
   hold on;
end