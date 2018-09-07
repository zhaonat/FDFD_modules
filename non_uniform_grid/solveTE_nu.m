
function [Hz, Ex, Ey, A,b, Dxf, Dyf, Dxb, Dyb] = ...
    solveTE_nu(L0, wvlen, xrange, yrange, eps_r, Mz, Npml,xs, ys)
%% Input Parameters
% wvlen: wavelength in L0
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% Mz: Nx-by-Ny array of magnetic current source density
% Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML
% xs, ys: 1D arrays containing the deltax and deltay gradings in the domain

%% Output Parameters
% Hz, Ex, Ey: Nx-by-Ny arrays of H- and E-field components

%% Set up the domain parameters.

%normal SI parameters
eps_0 = 8.85*10^-12*L0;
mu_0 = 4*pi*10^-7*L0; 
eps0 = eps_0;  % vacuum permittivity
mu0 = mu_0;  % vacuum permeability in
c0 = 1/sqrt(eps0*mu0);  % speed of light in 
N = size(eps_r);  % [Nx Ny] THIS IS THE POINT WHERE THE GRID SIZE IS DETERMINED
omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec

%% Set up the permittivity and permeability in the domain.
% bwdmean does nearest neighbor averaging (smoothes out stuff)
eps_x = bwdmean_w(eps0 * eps_r, 'y');  % average eps for eps_x
eps_y = bwdmean_w(eps0 * eps_r, 'x');  % average eps for eps_y
%these are fully dense matrices...

%% Set up number of cells
%the wavelength is much larger than the dimensions of the system...

Nx = N(1); dx = (xmax-xmin)/Nx; %we set dx and dy using the new xrange...
Ny = N(2); dy = (ymax-ymin)/Ny;
% Nz = 1; dz = 1; 2D solving only
M = prod(N); %total number of cells

%% Set up the Split coordinate PML
Nx_pml = Npml(1); Ny_pml = Npml(2);
sxf = create_sfactor_mine(xrange,'f',omega,eps_0,mu_0,Nx,Nx_pml);
syf = create_sfactor_mine(yrange,'f', omega,eps_0,mu_0,Ny,Ny_pml);
sxb = create_sfactor_mine(xrange, 'b', omega,eps_0,mu_0, Nx, Nx_pml);
syb = create_sfactor_mine(yrange,'b', omega,eps_0,mu_0,Ny,Ny_pml);

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
epxList = reshape(eps_x,M,1);
epyList = reshape(eps_y,M,1);
Tepx = spdiags(epxList,0,M,M); % creates an MxM matrix, which is the correct size,
%the M entries in epsList is put on the diagonals
Tepy = spdiags(epyList,0,M,M);
Tmz = mu0*speye(M); %in most cases, permeability is that of free-space

%% Create Magnetic vector Mz (source profile determined by Mz input)
Mz = reshape(Mz,M,1);
Mz = sparse(Mz);

%% create the derivative oeprators w/ PML

N = [Nx, Ny];
dL = [dx dy]; % Remember, everything must be in SI units beforehand

[Fsx, Fsy, Fsx_conj, Fsy_conj] = non_uniform_scaling(xs, ys);

Dxf = Fsx^-1*createDws('x', 'f', dL, N);%*Fsx; 
Dyf = Fsy^-1*createDws('y', 'f', dL, N);%*Fsy;
Dyb = Fsy_conj^-1*createDws('y', 'b', dL, N);%*Fsx_conj; 
Dxb = Fsx_conj^-1*createDws('x', 'b', dL, N);%*Fsy_conj; 
Dxf_pml = Sxf^-1*Dxf; 
Dyf_pml = Syf^-1*Dyf;
Dyb_pml = Syb^-1*Dyb; 
Dxb_pml = Sxb^-1*Dxb; 


%% Construct the matrix A, everything is in 2D
A = Dxf_pml*(Tepx^-1)*Dxb_pml + ...
    Dyf_pml*(Tepy^-1)*Dyb_pml + omega^2*Tmz;
%% construct the matrix b, everything is in 2D
b = 1i*omega*Mz;

%% solve system
% %% Solve the equation.
 if all(b==0)
 	hz = zeros(size(b));
 else
 	hz = A\b;
 end
 Hz = reshape(hz, N);

 %% now solve for Ex and Ey
 ey = (Tepx^-1*Dyb)*hz*1/(1i*omega);
 ex = 1/(1i*omega)*(Tepy^-1*Dxb)*hz;
 Ey = reshape(ey,N);
 Ex = reshape(ex,N);
 
end