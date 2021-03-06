

function [Hz, ky_eigs,A, Dxf, Dxb] = FDFD_1D_Ky_eigensolve(L0, dx, epsilon, omega, Nx, n, Npml)

% epsilon :  is an Nx by 1 array containing the dielectric profile of the waveguide array
% n = number of eigs

if(nargin < 7)
   Npml = 0;
end

%% parameter setup
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

N = [Nx, 1];
dL = [dx, 1];
Lx = Nx/dx;
xrange = [-Lx,Lx]/2;

Tep_x = spdiags(bwdmean_w(eps0*epsilon, 'x'),0,Nx,Nx); %still a Yee's grid kind of
Teps = spdiags(eps0*epsilon,0,Nx,Nx); 

%Tep = diag(eps0*epsilon);

%% get operators
Dxb = createDws('x', 'b', dL, N); 
Dxf = createDws('x', 'f', dL, N); 

%create pml
s_vector_x_f = create_sfactor(xrange, 'f', omega, eps0, mu0, Nx, Npml); 
s_vector_x_b = create_sfactor(xrange, 'b', omega, eps0, mu0, Nx, Npml); 
Sxf = diag(s_vector_x_f);
Sxb = diag(s_vector_x_b);

%create abc... only partially working
abc_f = ABC_1D(L0, Nx, dx, omega, 'f');
abc_b = ABC_1D(L0, Nx, dx, omega, 'b');

%add pml on
Dxb = Sxb\Dxb;
Dxf = Sxf\Dxf;
% Dxb = abc_b;
% Dxf = abc_f;

%% formulate equation
Derivative_op = Dxf*Tep_x^-1*Dxb;
Derivative_op(Nx,Nx-1) = 1; 
A = Teps*Derivative_op + Teps*omega^2*mu0;

% doing largestabs is not sufficient...
wvlen = 2*pi*c0/omega;
n_diel = sqrt(max(max(real(epsilon)))); 
beta_est = abs(2*pi*n_diel / wvlen); 


[vz_temp, ky_sqr] = eigs(A, n, beta_est^2); 
ky_eigs = sqrt(diag(ky_sqr));

for i = 1:n
    hz_temp = vz_temp(:, i);
%     ex_temp = 1/(1i*omega) * Tep^-1 * Dyb * hz_temp; 
%     ey_temp = 1/(1i*omega) * Tep^-1 * (-Dxb * hz_temp); 
    
    Hz{i} = hz_temp; 
%     Ex{i} = ex_temp;
%     Ey{i} = ey_temp;
end



end