function A = TEWKSolver(L0, xrange, yrange, eps_r, Npml,K, omega_guess)
    %normal SI parameters
    eps_0 = 8.85*10^-12*L0;
    mu_0 = 4*pi*10^-7*L0; 

    N = size(eps_r);  % [Nx Ny] THIS IS THE POINT WHERE THE GRID SIZE IS DETERMINED
    omega = omega_guess;
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
    epzList = reshape(eps_0*eps_r,M,1);
    Tepz = spdiags(epzList,0,M,M); % creates an MxM matrix, which is the correct size,
    Tmz = mu_0*speye(M); %in most cases, permeability is that of free-space



    %% create the derivative oeprators w/ PML

    N = [Nx, Ny];
    dL = [dx dy]; % Remember, everything must be in SI units beforehand
    Dxf = createDws('x', 'f', dL, N); 
    Dyf = createDws('y', 'f', dL, N);
    Dyb = createDws('y', 'b', dL, N); 
    Dxb = createDws('x', 'b', dL, N); 
    Dxf_pml = Sxf^-1*Dxf; 
    Dyf_pml = Syf^-1*Dyf;
    Dyb_pml = Syb^-1*Dyb; 
    Dxb_pml = Sxb^-1*Dxb; 

    Dxs_K = (Dxf_pml*Dxb_pml - 1i*K*Dxf_pml - 1i*K*Dxb_pml - K^2*speye(M));
    %% Construct the matrix A, everything is in 2D
    A = (Tmz*Tepz)^-1*(Dxs_K + Dyf_pml*Dyb_pml);
    % note a warning about ill-conditioned matrices will pop up here, but
    % for our purposes, it is okay.

 
end