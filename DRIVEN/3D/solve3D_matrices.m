function [A, b, Ao, bo, omega, TepsSuper, TmuSuper] = ...
    solve3D_matrices(L0, wvlen, xrange, yrange, zrange, eps_r, JCurrentVector, Npml, s)

    %% Set up the domain parameters.
    %% Ao, bo are the matrices without Wonseok's correction
    %% this function only produces the system, but does not solve anything
    
    %normal SI parameters
    eps_0 = 8.85*10^-12*L0;
    mu_0 = 4*pi*10^-7*L0; 
    eps0 = eps_0;  % vacuum permittivity
    mu0 = mu_0;  % vacuum permeability in
    c0 = 1/sqrt(eps_0*mu_0);  % speed of light in 
    N = size(eps_r);  % [Nx Ny] THIS IS THE POINT WHERE THE GRID SIZE IS DETERMINED
    omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec

    %% Set up the permittivity and permeability in the domain.
    % bwdmean does nearest neighbor averaging (smoothes out stuff)
    eps_x = bwdmean_w(eps0 * eps_r, 'y');  % average eps for eps_x
    eps_y = bwdmean_w(eps0 * eps_r, 'x');  % average eps for eps_y
    eps_z = bwdmean_w(eps0 * eps_r, 'z');
    %these are fully dense matrices...

    %currently, eps_x and eps_y are ultra-dense, which isn't right...

    %% Set up number of cells
    %the wavelength is much larger than the dimensions of the system...
    xmin = xrange(1); xmax = xrange(2);
    ymin = yrange(1); ymax = yrange(2);
    zmin = zrange(1); zmax = zrange(2);
    Nx = N(1); dx = (xmax - xmin)/Nx;
    Ny = N(2); dy = (ymax - ymin)/Ny;
    Nz = N(3); dz = (zmax - zmin)/Nz;
    % Nz = 1; dz = 1; 2D solving only
    M = prod([Nx, Ny, Nz]); %total number of cells
    Mx = JCurrentVector(1:Nx,:,:); My = JCurrentVector(Nx+1:2*Nx,:,:); Mz = JCurrentVector(2*Nx+1: 3*Nx,:,:);
    
    %% Create the dielectric and permeability arrays (ex, ey, muz)
    %create a diagonal block matrix of ep and mu...
    epxList = reshape(eps_x,M,1);
    epyList = reshape(eps_y,M,1);
    epzList = reshape(eps_z,M,1);
    Tepx = spdiags(epxList,0,M,M); % creates an MxM matrix, which is the correct size,
    %the M entries in epsList is put on the diagonals
    Tepy = spdiags(epyList,0,M,M);
    Tepz = spdiags(epzList, 0,M,M);
    Tmz = mu0*speye(M); %in most cases, permeability is that of free-space

    %% Create SuperVector of the Dielectrics and Permeability
    TepsSuper = blkdiag(Tepx,Tepy,Tepz);
    TmuSuper = blkdiag(Tmz, Tmz, Tmz);
    %% Create Current Source vector J 
    % dimension = M*1
    Mz = reshape(Mz, M, 1);
    My = reshape(My, M, 1);
    Mx = reshape(Mx, M, 1);
    Mz = sparse(Mz);

    %% create the derivative oeprators w/ PML
    Nx_pml = Npml(1); Ny_pml = Npml(2); Nz_pml = Npml(3);
    dL = [dx dy dz]; % Remember, everything must be in SI units beforehand

    %% PML operators
    sxf = create_sfactor_mine(xrange,'f',omega,eps_0,mu_0,Nx,Nx_pml);
    syf = create_sfactor_mine(yrange,'f', omega,eps_0,mu_0,Ny,Ny_pml);
    sxb = create_sfactor_mine(xrange, 'b', omega,eps_0,mu_0, Nx, Nx_pml);
    syb = create_sfactor_mine(yrange,'b', omega,eps_0,mu_0,Ny,Ny_pml);
    szf = create_sfactor_mine(xrange, 'f', omega,eps_0,mu_0, Nz, Npml(3));
    szb = create_sfactor_mine(yrange,'b', omega,eps_0,mu_0,Nz,Npml(3));


    % now we create the matrix (i.e. repeat sxf Ny times repeat Syf Nx times)
    [Sxf, Syf, Szf] = ndgrid(sxf, syf, szf);
    [Sxb, Syb, Szb] = ndgrid(sxb, syb, szb);

    %Sxf(:) converts from n x n t0 n^2 x 1
    Sxf=spdiags(Sxf(:),0,M,M);
    Sxb=spdiags(Sxb(:),0,M,M);
    Syf=spdiags(Syf(:),0,M,M);
    Syb=spdiags(Syb(:),0,M,M);
    Szf=spdiags(Szf(:),0,M,M);
    Szb=spdiags(Szb(:),0,M,M);
    
    
    %% Derivative operators
    Dxf = Sxf\createDws('x', 'f', dL, N); 
    Dyf = Syf\createDws('y', 'f', dL, N);
    Dyb = Syb\createDws('y', 'b', dL, N); 
    Dxb = Sxb\createDws('x', 'b', dL, N); 
    Dzf = Szf\createDws('z', 'f', dL, N); 
    Dzb = Szb\createDws('z', 'b', dL, N); 
    

    %% Construct Ch and Ce operators (curlH, curl E) or basically curl_forward
    zero_matrix = zeros(M,'like',Dxf) ;
    % and curl backward
    Ce = [zero_matrix -Dzf Dyf; Dzf zero_matrix -Dxf; -Dyf Dxf zero_matrix];
    Ch = [zero_matrix -Dzb Dyb; Dzb zero_matrix -Dxb; -Dyb Dxb zero_matrix];

    %% Construct the matrix A, everything is in 2D
    %% constrct the grad(eE) term
    %gradient(divergence)
    GradDiv = [Dxf*Tepx^-1*Dxb*Tepx, Dxf*Tepx^-1*Dyb*Tepy, Dxf*Tepx^-1*Dzb*Tepz; ...
        Dyf*Tepy^-1*Dxb*Tepx, Dyf*Tepy^-1*Dyb*Tepy, Dyf*Tepy^-1*Dzb*Tepz; ...
        Dzf*Tepz^-1*Dxb*Tepx, Dzf*Tepz^-1*Dyb*Tepy, Dzf*Tepz^-1*Dzb*Tepz];

    WAccelScal = speye(3*(Nx*Ny*Nz))*TmuSuper^-1;
    A = Ch*TmuSuper^-1*Ce + s*WAccelScal*GradDiv - omega^2*TepsSuper;
    Ao = Ch*TmuSuper^-1*Ce - omega^2*TepsSuper;

    %% construct the matrix b, everything is in 2D
    J = [Mx; My; Mz];
    b = -1i*omega*J; bo = b;
    JCorrection = (1i/omega) * (s*GradDiv*WAccelScal)*TepsSuper^-1*J;
    b = b+JCorrection;


end