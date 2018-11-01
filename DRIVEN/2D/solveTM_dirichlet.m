function [Ez, Hx, Hy, A, omega,b, Dxf, Dyf, Dxb, Dyb] = solveTM_dirichlet(L0, wvlen, xrange, yrange, eps_r, Jz, Npml)
%% Input Parameters
% L0: length unit (e.g., L0 = 1e-9 for nm)
% wvlen: wavelength in L0
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% Mz: Nx-by-Ny array of magnetic current source density
% Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML

%% Output Parameters
% Hz, Ex, Ey: Nx-by-Ny arrays of H- and E-field components
% dL: [dx dy] in L0
% A: system matrix of A x = b
% omega: angular frequency for given wvlen

%% Set up the domain parameters.

%normal SI parameters
eps0 = 8.85*10^-12*L0;
mu0 = 4*pi*10^-7*L0;
c0 = 1/sqrt(eps0*mu0);  % speed of light in 
N = size(eps_r);  % [Nx Ny]
L = [diff(xrange) diff(yrange)];  % [Lx Ly]
dL = L./N;  % [dx dy]
omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec


%% Set up number of cells
%the wavelength is much larger than the dimensions of the system...
xmin = xrange(1); xmax = xrange(2);
ymin = yrange(1); ymax = yrange(2);
Nx = N(1);
Ny = N(2); 
% Nz = 1; dz = 1; 2D solving only

M = prod([Nx, Ny]); %total number of cells

%% Set up the Split coordinate PML
%sx = create_sfactor('f',Nx);
%sy = creates_factor('f',Ny);
Nx_pml = Npml(1); Ny_pml = Npml(2);
Nwx = Nx; Nwy = Ny;
wrangex = [1, Nwx]; wrangey = [1, Nwy];
sxf = create_sfactor(wrangex,'f',omega,eps0,mu0,Nwx,Nx_pml);
syf = create_sfactor(wrangey,'f', omega,eps0,mu0,Nwy,Ny_pml);
sxb = create_sfactor(wrangex, 'b', omega,eps0, mu0, Nwx, Nx_pml);
syb = create_sfactor(wrangey,'b', omega,eps0,mu0,Nwy,Ny_pml);
% now we create the matrix
[Sxf, Syf] = ndgrid(sxf,syf);
[Sxb, Syb] = ndgrid(sxb, syb);


%Sxf(:) converts from n x n t0 n^2 x 1
Sxf=spdiags(Sxf(:),0,M,M);
Sxb=spdiags(Sxb(:),0,M,M);
Syf=spdiags(Syf(:),0,M,M);
Syb=spdiags(Syb(:),0,M,M);


%% Create the dielectric and permeability arrays (ex, ey, muz)
eps_z = bwdmean_w(eps0 *eps_r, 'z');

%create a diagonal block matrix of ep and mu...
epzList = reshape(eps_z,M,1);
Tepz = spdiags(epzList,0,M,M); % creates an MxM matrix, which is the correct size,
%the M entries in epsList is put on the diagonals
Tmz = mu0*speye(M); %in most cases, permeability is that of free-space
Tmy = Tmz; Tmx = Tmz;

%% Create Magnetic Current Source Mz at center of the grid
% dimension = M*1
Jz = reshape(Jz,M,1);
Jz = sparse(Jz);

%% create the derivative oeprators
Dxf = Sxf\createDws_dirichlet('x', 'f', dL, N); 
Dyf = Syf\createDws_dirichlet('y', 'f', dL, N);
Dyb = Syb\createDws_dirichlet('y', 'b', dL, N); 
Dxb = Sxb\createDws_dirichlet('x', 'b', dL, N);

%% Construct the matrix A, everything is in 2D
A = Dxb*(Tmy^-1)*Dxf + Dyb*(Tmx^-1)*Dyf + omega^2*Tepz;
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

 
end