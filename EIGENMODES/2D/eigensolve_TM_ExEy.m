

function [eigenvals, eigenmodes,A] = eigensolve_TM_ExEy(L0, wvlen, xrange, ...
    yrange, eps_tensor, Npml, neigs)
   
    %% Input Parameters
    % wvlen: wavelength in L0
    % xrange: [xmin xmax], range of domain in x-direction including PML
    % yrange: [ymin ymax], range of domain in y-direction including PML
    % eps_tensor: 2x2 cell array, each cell is an Nx x Ny array containing
    % the spatial grid with the epsilon tensor values
    % exx exy
    % eyx eyy
    % Mz: Nx-by-Ny array of magnetic current source density
    % Npml: [Nx_pml Ny_pml], number of cells in x- and y-normal PML
    % neigs: number of eigenmodes to get
    
    %% Output Parameters
    % Hz, Ex, Ey: Nx-by-Ny arrays of H- and E-field components
    % dL: [dx dy] in L0
    % A: system matrix of A x = b
    % omega: angular frequency for given wvlen

    %% Set up the domain parameters.

    %normal SI parameters
    eps_0 = 8.85*10^-12*L0;
    mu_0 = 4*pi*10^-7*L0; 
    eps0 = eps_0;  % vacuum permittivity
    mu0 = mu_0;  % vacuum permeability in
    c0 = 1/sqrt(eps0*mu0);  % speed of light in 
    omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec

    %% unravel the epsilon tensor;
    exx = eps_tensor{1,1};
    exy = eps_tensor{1,2};
    eyx = eps_tensor{2,1};
    eyy = eps_tensor{2,2};
    
    N = size(exx);
    
    %% Set up the permittivity and permeability in the domain.
    exx = bwdmean_w(eps0 * exx, 'y');  % average eps for eps_x
    eyy = bwdmean_w(eps0 * eyy, 'x');  % average eps for eps_y
    exy = bwdmean_w(eps0 * exy, 'y');  % average eps for eps_x
    eyx = bwdmean_w(eps0 * eyx, 'x');  % average eps for eps_y

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
    Exx = reshape(exx,M,1);
    Eyy = reshape(eyy,M,1);
    Exy = reshape(exy,M,1);
    Eyx = reshape(exy,M,1);
    
    Texx = spdiags(Exx,0,M,M); % creates an MxM matrix, which is the correct size,
    Teyy = spdiags(Eyy,0,M,M);
    Texy = spdiags(Exy, 0,M,M);
    Teyx = spdiags(Eyx, 0,M,M);
    
    Tep = eps0*[Texx, Texy; Teyx, Teyy];
    
    %% create the derivative oeprators w/ PML

    N = [Nx, Ny];
    dL = [dx dy]; % Remember, everything must be in SI units beforehand

    Dxf = createDws('x', 'f', dL, N); 
    Dyf = createDws_dense('y', 'f', dL, N);
    Dyb = createDws_dense('y', 'b', dL, N); 
    Dxb = createDws_dense('x', 'b', dL, N); 
    Dxf = Sxf^-1*Dxf; Dyf = Syf^-1*Dyf;
    Dyb = Syb^-1*Dyb; Dxb = Sxb^-1*Dxb; 
    Deriv = [-Dyb*Dyf Dyb*Dxf; Dxb*Dyf, -Dxb*Dxf];
    
    %% PROBLEM CONSTRUCTION
    disp('get operator');
    A = Tep\Deriv;
    
    disp('start eigensolve');
    %% eigensolver
    beta_est = 2*pi/(wvlen*L0);
    [U,V] = eigs(A, neigs, 'smallestabs');
    eigenvals = diag(V); %eigenvals solved are omega^2*mu0
    eigenvals = (1/mu0)*eigenvals;
    
    eigenmodes = cell(1);
    %% process the eigenmodes
    %V = Ex, Ey but flattened;
    for i = 1:neigs
        ex =U(1:M,i); ey = U( M+1:2*M,i);
        Ex = reshape(ex, Nx,Ny);
        Ey = reshape(ey, Nx, Ny);
        %reconstruct Hz
        Hz = (Dxf*ey - Dyf*ex)/(1i*omega);
        eigenmodes{i,1} = Ex;
        eigenmodes{i,2} = Ey;
        eigenmodes{i,3} = Hz;
        
    end

end