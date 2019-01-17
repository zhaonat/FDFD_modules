

function [Hz, ky_eigs,A, Dxf, Dxb, Tep_x, Teps] = ...
    FDFD_1D_Bloch_EH(L0, dx, epsilon, omega, Nx,Kx, n)

%% solves the eigenvalue problem simultaneously for E and H...
% allegedly hermitian for any real epsilon.

% epsilon :  is an Nx by 1 array containing the dielectric profile of the waveguide array
% n = number of eigs
%% THIS function solves a 1D system but it has a bloch boundary condition
% at the edges (note that the eigensolves does not solve for this 

%% parameter setup
eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec

N = [Nx, 1];
dL = [dx, 1];
Lx = Nx/dx;
xrange = [-Lx,Lx]/2;

Tep_x = spdiags( bwdmean_w(eps0*epsilon, 'x'),0,Nx,Nx); %still a Yee's grid kind of
Teps = spdiags(eps0*epsilon,0,Nx,Nx); 

%Tep = diag(eps0*epsilon);
L = [Lx,0]; %only creating Dxf and Dxb so we can having 0 in the position of Ly
k_vec = [Kx,0]; %only implement an x_bloch boundary...
%% get operators
Dxb = createDws_bloch('x', 'b', dL, N, k_vec, L); 
Dxf = createDws_bloch('x', 'f', dL, N, k_vec, L); 

I = speye(M);
%% formulate equation
A1 = (c0/(omega*mu0))*Dyf*Dxb + Teps*omega*mu0;
A2 = I*mu0*omega/c0 + (c0/omega)*Dyf*Teps_x^-1*Dyb;
A = blkdiag(A1,A2);
B = blkdiag()

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