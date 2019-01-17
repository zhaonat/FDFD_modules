function [Hz_modes, Ex_modes, Ey_modes, eigenvals,C,B,A, LHS, RHS] = ...
    eigensolve_TM_dispersive_anisotropic_Kx(L0, wvlen,...
    xrange, yrange, eps_tensor, Npml,neigs, kx_guess)
%% Input Parameters
% wvlen: wavelength in L0
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_tensor: 2x2 anisotropic (potentially permittivity)
% Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML

%% Output Parameters

%% Set up the domain parameters.

%normal SI parameters
eps_0 = 8.85*10^-12*L0;
mu_0 = 4*pi*10^-7*L0; 
eps0 = eps_0;  % vacuum permittivity
mu0 = mu_0;  % vacuum permeability in
c0 = 1/sqrt(eps0*mu0);  % speed of light in 
N = size(eps_tensor{1,1});  % [Nx Ny] THIS IS THE POINT WHERE THE GRID SIZE IS DETERMINED
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

%% create the derivative oeprators w/ PML

N = [Nx, Ny];
dL = [dx dy]; % Remember, everything must be in SI units beforehand

Dxf = createDws('x', 'f', dL, N); 
Dyf = createDws('y', 'f', dL, N);
Dyb = createDws('y', 'b', dL, N); 
Dxb = createDws('x', 'b', dL, N); 
Dxf = Sxf\Dxf; 
Dyf = Syf\Dyf;
Dyb = Syb\Dyb; 
Dxb = Sxb\Dxb; 

%% construct the inverse epsilon terms
fxx = TD*Texx;
fyy = TD*Teyy;
fxy = TD*Texy;
fyx = TD*Teyx;

%% Construct the matrix A, everything is in 2D
%% rule of thumb b is before f since we're operating on H
I = speye(M);
Z = zeros(M,M,'like',I);
C = -mu0^-1*fxx;                                                     %% K^2
B = mu0^-1*1i*(fxx*Dxb + Dxf*fxx + Dyf*fxy + fxy *Dyb);                           %% K
A = mu0^-1*(Dxf*fxx*Dxb  + Dxf*fyx*Dyb +...
    + Dyf*fxy*Dxb + Dyf*fyy*Dyb + omega^2*mu0*I);      %% CONSTANT

%% formulate the problem
LHS = [Z , -C; 
       I ,  Z];
RHS = [A, B; 
       Z, I ];
NL = size(LHS);

%% solve system
tic
[U,V] = eigs(RHS,LHS,neigs,kx_guess); %polyeigs is actually not suitable because it will solve all eigens
toc
%we need to linearize the eigenproblem ourselves

%% EXTRACT MODES
eigenvals = diag(V); %eigenvals solved are omega^2*mu0
Hz_modes = cell(1);
Ex_modes = cell(1);
Ey_modes = cell(1);
uz_bloch_modes = cell(1);
ux_bloch_modes = cell(1);
uy_bloch_modes = cell(1);

%% process the eigenmodes

x = linspace(xrange(1), xrange(2), N(1))+diff(xrange)/2; 
%the fact that xrange is positive and negative in x will be a problem
% when complex wavevectors are involved.

x = repmat(x.', [N(2), 1]); 
for i = 1:neigs
    Kx = V(i,i);
    hz = U(1:round(NL/2),i);%.*exp(-1i*real(Kx)*x);%.*exp(-1i*Kx*x);
    Hz = reshape(hz, Nx,Ny);
    uz_bloch_modes{i} = reshape(U(1:round(NL/2),i), Nx,Ny);
    
    ex = (1/(1i*omega))*TD*(Teyy*Dyb+Texy*Dxb)*hz;
    ey = -(1/(1i*omega))*TD*(Teyx*Dyb+Texx*Dxb)*hz;
    
    Ex = reshape(ex, Nx, Ny);
    Ey = reshape(ey, Nx,Ny);
    Hz_modes{i} = Hz;
    Ex_modes{i} = Ex;
    Ey_modes{i} = Ey;

end


 
end