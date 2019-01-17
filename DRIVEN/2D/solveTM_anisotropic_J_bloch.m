function [Hz, Ex, Ey, A, A_mode, omega,b] = ...
    solveTM_anisotropic_J_bloch(L0, wvlen, xrange, yrange, eps_tensor, Mz,Jx,Jy, Npml)
%% Input Parameters
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
% Dxf, Dyf - derivative matrices
% SXf, sxf - derivative matrices with pml implemented

%% Set up the domain parameters.

%normal SI parameters
eps_0 = 8.85*10^-12*L0;
mu_0 = 4*pi*10^-7*L0; 
eps0 = eps_0;  % vacuum permittivity
mu0 = mu_0;  % vacuum permeability in
c0 = 1/sqrt(eps0*mu0);  % speed of light in 
N = size(Mz);  % [Nx Ny] THIS IS THE POINT WHERE THE GRID SIZE IS DETERMINED
omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec



%% Set up number of cells
%the wavelength is much larger than the dimensions of the system...
xmin = xrange(1); xmax = xrange(2);
ymin = yrange(1); ymax = yrange(2);
Nx = N(1); dx = (xmax-xmin)/Nx;
Ny = N(2); dy = (ymax-ymin)/Ny;
% Nz = 1; dz = 1; 2D solving only
M = prod(N); %total number of cells

%% Set up the Split coordinate PML
%sx = create_sfactor('f',Nx);
%sy = creates_factor('f',Ny);
Nx_pml = Npml(1); Ny_pml = Npml(2);
sxf = create_sfactor_mine(xrange,'f',omega,eps_0,mu_0,Nx,Nx_pml, -12, 3.5);
syf = create_sfactor_mine(yrange,'f', omega,eps_0,mu_0,Ny,Ny_pml,-12, 3.5);
sxb = create_sfactor_mine(xrange, 'b', omega,eps_0,mu_0, Nx, Nx_pml,-12, 3.5);
syb = create_sfactor_mine(yrange,'b', omega,eps_0,mu_0,Ny,Ny_pml,-12, 3.5);

% now we create the matrix (i.e. repeat sxf Ny times repeat Syf Nx times)
[Sxf, Syf] = ndgrid(sxf, syf);
[Sxb, Syb] = ndgrid(sxb, syb);


%Sxf(:) converts from n x n t0 n^2 x 1
Sxf=spdiags(Sxf(:),0,M,M);
Sxb=spdiags(Sxb(:),0,M,M);
Syf=spdiags(Syf(:),0,M,M);
Syb=spdiags(Syb(:),0,M,M);


%% Create the dielectric and permeability arrays (ex, ey, muz)
exx = eps0*eps_tensor{1,1}; % this has not yet been flattened out into an NM by NM matrix
exy = eps0*eps_tensor{1,2};
eyx = eps0*eps_tensor{2,1};
eyy = eps0*eps_tensor{2,2};

%% Set up the permittivity and permeability in the domain.
% this should be x and y too for exx and eyy respectively 
exx = bwdmean_w(exx, 'y');  % average eps for eps_x
eyy = bwdmean_w(eyy, 'x');  % average eps for eps_y
exy =  bwdmean_w(exy, 'y');  % average eps for eps_x
eyx = bwdmean_w(eyx, 'x');  % average eps for eps_y
determinant = (exx.*eyy -exy.*eyx).^-1;

% flatten out
Exx = reshape(exx,M,1);
Eyy = reshape(eyy,M,1);
Exy = reshape(exy,M,1);
Eyx = reshape(eyx,M,1);
Determinant = reshape(determinant, M,1);

%% create interpolation functions;
Rxf = interpolate(N, 'x', 'f');
Rxb = interpolate(N, 'x', 'b');
Ryf = interpolate(N, 'y', 'f');
Ryb = interpolate(N, 'y', 'b');

% convert to matrices
Texx = spdiags(Exx,0,M,M); % creates an MxM matrix, which is the correct size,
Teyy = spdiags(Eyy,0,M,M);
%have to be careful...only invert nonzero exy or eyx
Texy = Ryb*Rxf*spdiags(Exy, 0,M,M);
Teyx = Rxb*Ryf*spdiags(Eyx, 0,M,M); %have large regions of the diagonal that are 0
TD = spdiags(Determinant, 0,M,M);

Tmz = mu0*speye(M); %in most cases, permeability is that of free-space

%% Create Magnetic vector Mz (source profile determined by Mz input)
% dimension = M*1

Mz = reshape(Mz,M,1);
Mz = sparse(Mz);


%% create the derivative oeprators w/ PML

N = [Nx, Ny];
dL = [dx dy]; % Remember, everything must be in SI units beforehand

Dxf = createDws('x', 'f', dL, N); 
Dyf = createDws('y', 'f', dL, N);
Dyb = createDws('y', 'b', dL, N); 
Dxb = createDws('x', 'b', dL, N); 
Dxf_pml = Sxf\Dxf; 
Dyf_pml = Syf\Dyf;
Dyb_pml = Syb\Dyb; 
Dxb_pml = Sxb\Dxb; 

%% Work on the electric current source
Jxf = sparse(reshape(Jx, M,1));
Jyf = sparse(reshape(Jy, M,1));

jxcon = Dyf_pml*(TD*(Teyy*Jxf+Texy*Jyf));
jycon = Dxf_pml*(TD*(Teyx*Jxf+Texx*Jyf));

%% Construct the matrix A, everything is in 2D
A_mode = Dyf_pml*TD*(Teyy*Dyb_pml + Texy *Dxb_pml) + ...
         Dxf_pml*TD*(Teyx*Dyb_pml + Texx*Dxb_pml);
A = A_mode + omega^2*Tmz;
% note a warning about ill-conditioned matrices might pop up here, but
% for our purposes, it is okay.

%% construct the matrix b, everything is in 2D
b = 1i*omega*Mz +jxcon +jycon;

%% solve system
% %% Solve the equation.
if all(b==0)
hz = zeros(size(b));
else
hz = A\b;
end
Hz = reshape(hz, N);

%% reconstruct Ex Ey
ex = TD*(Teyy*Dyb_pml+Texy*Dxb_pml)*hz;
ey = -TD*(Teyx*Dyb_pml+Texx*Dxb_pml)*hz;
Ex = (1/(1i*omega))*reshape(ex, N);
Ey = (1/(1i*omega))*reshape(ey, N);

 
end