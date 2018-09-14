function [Hz, Ex, Ey, A, A_mode, b, Tepx, Tepy, Dyb, Dxb] = ...
    solveTE_Jdipole(L0, wvlen, xrange, yrange, eps_r, Mz, Npml, Jx, Jy)
%% Input Parameters
% L0: length unit (e.g., L0 = 1e-9 for nm)
% wvlen: wavelength in L0
% xrange: [xmin xmax], range of domain in x-direction including PML
% yrange: [ymin ymax], range of domain in y-direction including PML
% eps_r: Nx-by-Ny array of relative permittivity
% Mz: Nx-by-Ny array of magnetic current source density
% Jx, Jy: Nx-by-Ny array of ELECTRIC current source
% Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML

%% Output Parameters
% Ez, Hx, Hy: Nx-by-Ny arrays of H- and E-field components
% dL: [dx dy] in L0
% A: system matrix of A x = b
% omega: angular frequency for given wvlen

%% Set up the domain parameters.
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

N = size(eps_r);  % [Nx Ny]
L = [diff(xrange) diff(yrange)];  % [Lx Ly]
dL = L./N;  % [dx dy]

M = prod(N); 
Nx = N(1); Ny = N(2);
omega = 2*pi*c0/wvlen;  % angular frequency in rad/sec

%% Deal with the s_factor
Nx_pml = Npml(1); Ny_pml = Npml(2);
sxf = create_sfactor(xrange,'f',omega,eps0,mu0,Nx,Nx_pml);
syf = create_sfactor(yrange,'f', omega,eps0,mu0,Ny,Ny_pml);
sxb = create_sfactor(xrange, 'b', omega,eps0,mu0, Nx, Nx_pml);
syb = create_sfactor(yrange,'b', omega,eps0,mu0,Ny,Ny_pml);

% now we create the matrix (i.e. repeat sxf Ny times repeat Syf Nx times)
[Sxf, Syf] = ndgrid(sxf, syf);
[Sxb, Syb] = ndgrid(sxb, syb);


%Sxf(:) converts from n x n t0 n^2 x 1
Sxf=spdiags(Sxf(:),0,M,M);
Sxb=spdiags(Sxb(:),0,M,M);
Syf=spdiags(Syf(:),0,M,M);
Syb=spdiags(Syb(:),0,M,M);

%% Set up the permittivity and permeability in the domain.
eps_x = bwdmean_w(eps0*eps_r, 'x'); 
eps_y = bwdmean_w(eps0*eps_r, 'y'); 

% Setup the Teps_x, Teps_y, and Tmu_z matrices
Tepx = spdiags(eps_x(:), 0, M, M); 
Tepy = spdiags(eps_y(:), 0, M, M); 


%% Construct derivate matrices
Dyb = Syb\createDws('y', 'b', dL, N); 
Dxb = Sxb\createDws('x', 'b', dL, N); 
Dxf = Sxf\createDws('x', 'f', dL, N); 
Dyf = Syf\createDws('y', 'f', dL, N); 

%% Reshape Mz into a vector
mz = reshape(Mz, M, 1); 
jx = reshape(Jx, M, 1);
jy = reshape(Jy, M, 1);

%% contribution of j
jx_con = -Dyf*(Tepx\jx);
jy_con = Dxf*(Tepy\jy);

%% Construct A matrix and b vector
% A = Sx_f *Dxf* T_eps_y^-1 *Sx_b *Dxb + Sy_f *Dyf* T_eps_x^-1* Sy_b*Dyb + omega^2*T_mu_z; 
% A = Sx_f*Dxf* T_eps_y^-1 *Sx_b*Dxb + Sy_f*Dyf* T_eps_x^-1* Sy_b*Dyb + omega^2*T_mu_z; 

I = speye(M); 
A_mode = Dxf*(Tepx^-1)*Dxb + Dyf*(Tepy^-1)*Dyb;
A = A_mode + omega^2*mu0*I; 
% % A = Dxf* T_eps_y^-1 *Dxb + Dyf* T_eps_x^-1* Dyb + omega^2*T_mu_z; 
b = 1i * omega * mz + jx_con + jy_con; 



%% Solve the equation.
if all(b==0)
	hz = zeros(size(b));
else
	hz = A\b;
end
Hz = reshape(hz, N);


ex = 1/(1i*omega) * Tepy^-1 * Dyb * hz; 
ey = 1/(1i*omega) * Tepx^-1 * (-Dxb * hz); 

Ex = reshape(ex, N); 
Ey = reshape(ey, N); 

end
