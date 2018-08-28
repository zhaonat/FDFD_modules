
function [A, omega,b, Dzf, Dzb, Tep] = ...
    solveFDFD_1D(wvlen, zrange, dL, eps_r, Mz, Npml, L0)

    eps_0 = 8.85*10^-12*L0;
    mu_0 = 4*pi*10^-7*L0; 
    eps0 = eps_0;  % vacuum permittivity
    mu0 = mu_0;  % vacuum permeability in
    c0 = 1/sqrt(eps0*mu0);  % speed of light in 
    N = [length(eps_r), 1];  
    omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec

    %% Set up number of cells
    %the wavelength is much larger than the dimensions of the system...
    zmin = zrange(1); zmax = zrange(2);
    Nx = N(1); 
    % Nz = 1; dz = 1; 2D solving only
    
    %% Generate PML in z direction
    
    %% reshape dielectric
    Tep = eps0*diag(eps_r);
    
    %% create Dzf
    Dzf = createDws_dense('x', 'f', dL, N); 
    Dzb = createDws_dense('x', 'b', dL, N);
    
    %% formulate A operator
    A = Dzf*Tep^-1*Dzb + omega^2*mu0*eye(Nx);
    
    %% construct the source b
    b = 1i*omega*Mz;
 
end
