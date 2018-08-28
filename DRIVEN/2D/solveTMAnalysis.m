%% SOLVE TM ANALYSIS
%% Parameter Setup
c0 = 3*10^8;
L0 = 1e-6;  % length unit: microns
wvlen = 4.0*L0;  % wavelength in L0
xrange = [-5 5]*L0;  % x boundaries in L0
yrange = [-5 5]*L0;  % y boundaries in L0
Nx = 100; Ny = 100;
N = [Nx Ny];  % [Nx Ny]
Npml = [10 10];  % [Nx_pml Ny_pml]
mu0 = 4*pi*10^-7; mu_0 = mu0; mu = mu0;
eps0 = 8.85*10^-12; eps_0 = eps0;
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
eps_r = ones(N);

%% Fancy Structure

% eps_r(Nx/2:Nx/2+10,:) = 12;
% for j = 1:Ny
%     if(mod(j,10) == 0 || mod(j,11) == 0 || mod(j, 12) == 0)
%        eps_r(Nx/2:Nx/2+10,j) = 0; 
%     end
% end

%% Set up the 1agnetic current source density.
Jz = zeros(N);
ind_src = [30 30];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Jz(ind_src(1), ind_src(2)) = 1;


%% SOLVER CODE STARTS HERE
%normal SI parameters
eps_0 = 8.85*10^-12;
mu_0 = 4*pi*10^-7; 
eps0 = eps_0;  % vacuum permittivity
mu0 = mu_0;  % vacuum permeability in
c0 = 1/sqrt(eps0*mu0);  % speed of light in 
N = size(eps_r);  % [Nx Ny] THIS IS THE POINT WHERE THE GRID SIZE IS DETERMINED
omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec

%% Set up the permittivity and permeability in the domain.
% bwdmean does nearest neighbor averaging (smoothes out stuff)

eps_z = bwdmean_w(eps0 *eps_r, 'z');
%these are fully dense matrices...

%currently, eps_x and eps_y are ultra-dense, which isn't right...

%% Set up number of cells
%the wavelength is much larger than the dimensions of the system...
xmin = xrange(1); xmax = xrange(2);
ymin = yrange(1); ymax = yrange(2);
Nx = N(1); dx = (xmax-xmin)/Nx;
Ny = N(2); dy = (ymax-ymin)/Ny;
% Nz = 1; dz = 1; 2D solving only
M = prod([Nx, Ny]); %total number of cells

%% Set up the Split coordinate PML
%sx = create_sfactor('f',Nx);
%sy = creates_factor('f',Ny);
Nx_pml = Npml(1); Ny_pml = Npml(2);
Nwx = Nx; Nwy = Ny;
sxf = create_sfactor_mine(xrange,'f',omega,eps_0,mu_0,Nwx,Nx_pml);
syf = create_sfactor_mine(yrange,'f', omega,eps_0,mu_0,Nwy,Ny_pml);
sxb = create_sfactor_mine(xrange, 'b', omega,eps_0,mu_0, Nwx, Nx_pml);
syb = create_sfactor_mine(yrange,'b', omega,eps_0,mu_0,Nwy,Ny_pml);

% now we create the matrix (i.e. repeat sxf Ny times repeat Syf Nx times)
[Sxf, Syf] = ndgrid(sxf, syf);
[Sxb, Syb] = ndgrid(sxb, syb);

%Sxf(:) converts from n x n t0 n^2 x 1
Sxf=spdiags(Sxf(:),0,M,M);
Sxb=spdiags(Sxb(:),0,M,M);
Syf=spdiags(Syf(:),0,M,M);
Syb=spdiags(Syb(:),0,M,M);


%% Create the dielectric and permeability arrays (ex, ey, muz)
%create a diagonal block matrix of ep and mu...
epzList = reshape(eps_z,M,1);
Tepz = spdiags(epzList,0,M,M); % creates an MxM matrix, which is the correct size,
%the M entries in epsList is put on the diagonals
Tmz = mu0*speye(M); %in most cases, permeability is that of free-space
Tmy = Tmz; Tmx = Tmz;

%% Create Magnetic vector Mz (source profile determined by Mz input)
% dimension = M*1
Jz = reshape(Jz,M,1);
Jz = sparse(Jz);

%% create the derivative oeprators w/ PML

N = [Nx, Ny];
dL = [dx dy]; % Remember, everything must be in SI units beforehand

Dxf = createDws_dense('x', 'f', dL, N); 
Dyf = createDws_dense('y', 'f', dL, N);
Dyb = createDws_dense('y', 'b', dL, N); 
Dxb = createDws_dense('x', 'b', dL, N); 
Dxf_pml = Sxf^-1*Dxf; 
Dyf_pml = Syf^-1*Dyf;
Dyb_pml = Syb^-1*Dyb; 
Dxb_pml = Sxb^-1*Dxb; 


%% Construct the matrix A, everything is in 2D
A = Dxb_pml*(Tmy^-1)*Dxf_pml + Dyb_pml*(Tmx^-1)*Dyf_pml + omega^2*Tepz;
% note a warning about ill-conditioned matrices will pop up here, but
% for our purposes, it is okay.

%% construct the matrix b, everything is in 2D
b = 1i*omega*Jz;

%% solve system
t0 =cputime
% %% Solve the equation.
 if all(b==0)
 	ez = zeros(size(b));
 else
   %hz = A\b;
 	ez = A\b;
 end
 trun = cputime-t0;
 Ez = reshape(ez, N);

 %% now solve for Ex and Ey
 hx = -1/(1i*omega)*(Tmx^-1*Dyf)*ez;
 hy = (Tmy^-1*Dxf)*ez*(1/(1i*omega));
 Hy = reshape(hy,N);
 Hx = reshape(hx,N);
 
 
 %% Function Test
[Ez, Hx, Hy, A, omega,b, Sxf, Dxf, Dyf, sxf, syf,trun] = solveTM(wvlen, xrange, yrange, eps_r, Jz, Npml)
 
 figure;
 visabs((Ez), xrange, yrange)
 figure;
 visabs((Hy), xrange, yrange)
 figure;
 visabs((Hx), xrange, yrange)
 
 moviereal(Ez, xrange, yrange)