
function [Hz_modes, Ex_modes, Ey_modes, eigenvals] = eigensolve_TM(L0, wvlen, xrange, ...
    yrange, eps_iso, Npml, neigs)
   
    %EIGENSOLVE_TM Summary of this function goes here
    %   Detailed explanation goes here
    eps_0 = 8.85*10^-12*L0;
    mu_0 = 4*pi*10^-7*L0; 
    eps0 = eps_0;  % vacuum permittivity
    mu0 = mu_0;  % vacuum permeability in
    c0 = 1/sqrt(eps0*mu0);  % speed of light in 
    omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec

    %% unravel the epsilon tensor;
    N = size(eps_iso);
    
    %% Set up the permittivity and permeability in the domain.
    exx = bwdmean_w(eps_0*eps_iso, 'y');  % average eps for eps_x
    eyy = bwdmean_w(eps_0*eps_iso, 'x');  % average eps for eps_y

    %% Set up number of cells
    xmin = xrange(1); xmax = xrange(2);
    ymin = yrange(1); ymax = yrange(2);
    Nx = N(1); dx = (xmax-xmin)/Nx;
    Ny = N(2); dy = (ymax-ymin)/Ny;
    M = prod([Nx, Ny]); %total number of cells

    %% Set up the Split coordinate PML
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
    Tex = spdiags(reshape(exx,M,1),0,M,M);
    Tey = spdiags(reshape(eyy, M,1),0,M,M);
    

    %% create the derivative oeprators w/ PML
    N = [Nx, Ny];
    dL = [dx dy]; % Remember, everything must be in SI units beforehand

    Dxf = createDws('x', 'f', dL, N); 
    Dyf = createDws('y', 'f', dL, N);
    Dyb = createDws('y', 'b', dL, N); 
    Dxb = createDws('x', 'b', dL, N); 
    Dxf = Sxf^-1*Dxf; Dyf = Syf^-1*Dyf;
    Dyb = Syb^-1*Dyb; Dxb = Sxb^-1*Dxb; 
    
    %% PROBLEM CONSTRUCTION
    disp('get operator');
    A = -(Dxf*Tex^-1*Dxb + Dyf*Tey^-1*Dyb); %
        
    %% eigensolver
    disp('start eigensolve');
    omega_est = 2*pi*c0/wvlen;
    %[U,V] = eigs(A, neigs, 'smallestabs');
    %find eigenmodes near desired frequency
    [U,V] = eigs(A, neigs, omega_est^2*mu0);

    eigenvals = diag(V); %eigenvals solved are omega^2*mu0
    
    Hz_modes = cell(1);
    Ex_modes = cell(1);
    Ey_modes = cell(1);
    %% process the eigenmodes
    for i = 1:neigs
        hz = U(:,i);
        omega = V(i,i);
        Hz = reshape(hz, Nx,Ny);
        ex = (1/(1i*omega))*Tex^-1*Dyb*hz;
        ey = -(1/(1i*omega))*Tey^-1*Dxb*hz;
        Ex = reshape(ex, Nx, Ny);
        Ey = reshape(ey, Nx,Ny);
        Hz_modes{i} = Hz;
        Ex_modes{i} = Ex;
        Ey_modes{i} = Ey;

    end
    
    
end

