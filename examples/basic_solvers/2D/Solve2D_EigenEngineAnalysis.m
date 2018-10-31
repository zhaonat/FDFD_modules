%% this explicitly sets the Nz dimension to be 0.
L0 = 1e-6;  % length unit: microns
wvlen = 4.0;  % wavelength in L0
xrange = [-5 5];  % x boundaries in L0
yrange = [-5 5];  % y boundaries in L0
Nx = 100; Ny = 100;
N = [Nx Ny];  % [Nx Ny]
M = N(1)*N(2);
Npml = [0 0];  % [Nx_pml Ny_pml]

%[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

Mz = zeros(N); My = Mz; Mx = Mz;
ind_src = [Nx/2 Ny/2];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;
JCurrentVector = [Mx; My; Mz];

s = 1;
%% Set up the permittivity.
eps_r = ones(N);

%% SOLVER CODE BEGINS HERE
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

    %% Set up the domain parameters.

    %normal SI parameters
    eps_0 = 8.85*10^-12;
    mu_0 = 4*pi*10^-7; 
    eps0 = eps_0;  % vacuum permittivity
    mu0 = mu_0;  % vacuum permeability in
    c0 = 1/sqrt(eps0*mu0);  % speed of light in 
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
    Nx = N(1); dx = (xmax-xmin)/Nx;
    Ny = N(2); dy = (ymax-ymin)/Ny;
    % Nz = 1; dz = 1; 2D solving only
    M = prod([Nx, Ny]); %total number of cells

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
    epxList = reshape(eps_x,M,1);
    epyList = reshape(eps_y,M,1);
    epzList = reshape(eps_z,M,1);
    Tepx = spdiags(epxList,0,M,M); % creates an MxM matrix, which is the correct size,
    %the M entries in epsList is put on the diagonals
    Tepy = spdiags(epyList,0,M,M);
    Tepz = spdiags(epzList, 0,M,M);
    Tmz = mu0*speye(M); %in most cases, permeability is that of free-space
    TepsSuper = blkdiag(Tepx,Tepy,Tepz);
    TmuSuper = blkdiag(Tmz, Tmz, Tmz);
    %% Create Magnetic vector Mz (source profile determined by Mz input)
    % dimension = M*1
    Mx = JCurrentVector(1:Nx,:,:); My = JCurrentVector(Nx+1:2*Nx,:,:); Mz = JCurrentVector(2*Nx+1: 3*Nx,:,:);
    Mz = reshape(Mz, M, 1);
    My = reshape(My, M, 1);
    Mx = reshape(Mx, M, 1);
    Mz = sparse(Mz);
    J = [Mx; My; Mz];
   
    %% create the derivative oeprators w/ PML

    N = [Nx, Ny];
    dL = [dx dy]; % Remember, everything must be in SI units beforehand

    Dxf = createDws_dense('x', 'f', dL, N); 
    Dyf = createDws_dense('y', 'f', dL, N);
    Dyb = createDws_dense('y', 'b', dL, N); 
    Dxb = createDws_dense('x', 'b', dL, N); 
    Dxf_pml = Sxf^-1*Dxf; 
    Dyf_pml = Syf^-1*Dyf;
    Dyb_pml = Syb^-1*Dyb; 
    Dxb_pml = Sxb^-1*Dxb; 
    
    %% Construct Wonsoek's Accelerator
    GradDiv = [Dxf_pml*Dxf_pml Dxf_pml*Dyf_pml sparse(M,M); Dyf_pml*Dxf_pml Dyf_pml*Dyf_pml sparse(M,M); ...
        sparse(M,M), sparse(M,M), sparse(M,M)];
    
    
    Ce = [sparse(M,M) sparse(M,M) Dyf_pml; sparse(M,M) sparse(M,M) -Dxf_pml; -Dyf_pml Dxf_pml sparse(M,M)];
    Ch = [sparse(M,M) sparse(M,M) Dyb_pml; sparse(M,M) sparse(M,M) -Dxb_pml; -Dyb_pml Dxb_pml sparse(M,M)];
    
    WAccelScal = speye(3*M) * TmuSuper^-1;

    %% PROBLEM CONSTRUCTION
    b = -1i*omega*J;
    JCorrection = s*(1i/omega) * GradDiv*WAccelScal*J;
    b = b+JCorrection;
    A = Ch*TmuSuper^-1*Ce - s*WAccelScal*GradDiv - omega^2*TepsSuper;
    solution = A\b;
    
    %% SOLUTION EXTRAction
    solLength = length(solution);
    
    Ex = solution(1:solLength/3);
    Ey = solution(solLength/3+1:solLength*(2/3));
    Ez = solution(solLength*(2/3)+1: solLength);
    h = (1i/omega)*TmuSuper^-1*(Ce*solution);
    
    Hx = h(1:solLength/3);
    Hy = h(solLength/3+1:solLength*(2/3));
    Hz = h(solLength*(2/3)+1: solLength);
    
    
 %% SOLVER CODE ENDS
 %% BEGIN VISUALIZING
Exs = reshape(Ex, Nx, Ny);
Eys = reshape(Ey, Nx, Ny);
Ezs = reshape(Ez, Nx, Ny);
Hxs = reshape(Hx, Nx, Ny);
Hys = reshape(Hy, Nx, Ny);
Hzs = reshape(Hz, Nx, Ny);

imagesc(abs(Ezs)) %% the other two E-fields should be 0.
