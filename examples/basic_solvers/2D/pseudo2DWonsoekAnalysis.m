clear all
close all

%% this sets Nz = 1;

L0 = 1e-6;  % length unit: microns
wvlen = 5.0*L0;  % wavelength in L0
xrange = [-5 5]*L0;  % x boundaries in L0
yrange = [-5 5]*L0;  % y boundaries in L0
zrange = [-0.001 0.001]*L0;
Nx = 40; Ny = 40; Nz = 2;
N = [Nx Ny Nz];  % [Nx Ny]
M = N(1)*N(2)*N(3);
Npml = [0 0 0];  % [Nx_pml Ny_pml]

%[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

Mz = zeros(N); My = Mz; Mx = Mz;
ind_src = [5 1 1];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2), ind_src(3)) = 1;
JCurrentVector = [Mx; My; Mz];


%% Set up the permittivity.
eps_r = 1*ones(N);
%% add structure
% eps_r(7:13, 7:13, 7:13) = 12;
%% SOLVER CODE BEGINS HERE 

    %% Set up the domain parameters.

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

    N = [Nx, Ny, Nz];
    dL = [dx dy dz]; % Remember, everything must be in SI units beforehand

    Dxf = createDws_dense('x', 'f', dL, N); 
    Dyf = createDws_dense('y', 'f', dL, N);
    Dyb = createDws_dense('y', 'b', dL, N); 
    Dxb = createDws_dense('x', 'b', dL, N); 
    Dzf = createDws_dense('z', 'f', dL, N); 
    Dzb = createDws_dense('z', 'b', dL, N); 

    %% Construct Ch and Ce operators
    Ce = [zeros(M) -Dzf Dyf; Dzf zeros(M) -Dxf; -Dyf Dxf zeros(M)];
    Ch = [zeros(M) -Dzb Dyb; Dzb zeros(M) -Dxb; -Dyb Dxb zeros(M)];

    %% Construct the matrix A, everything is in 2D
    s = -1;
    %% constrct the eE term
    %gradient(divergence)
    Delf=[Dxf, Dyf, Dzf]; Delb = [Dxb, Dyb, Dzb];
    GradDiv = transpose(Delf)*Delb;
    
    GradDivHom = [Dxf*Dxb, Dxf*Dyb, Dxf*Dzb; ...
        Dyf*Dxb, Dyf*Dyb, Dyf*Dzb; ...
        Dzf*Dxb, Dzf*Dyb, Dzf*Dzb];
    
    GradDivMan = [Dxf*Tepx^-1*Dxb*Tepx, Dxf*Tepx^-1*Dyb*Tepy, Dxf*Tepx^-1*Dzb*Tepz; ...
        Dyf*Tepy^-1*Dxb*Tepx, Dyf*Tepy^-1*Dyb*Tepy, Dyf*Tepy^-1*Dzb*Tepz; ...
        Dzf*Tepz^-1*Dxb*Tepx, Dzf*Tepz^-1*Dyb*Tepy, Dzf*Tepz^-1*Dzb*Tepz];
    
    GradDivCurrent = [Dxf*Tepx^-1*Dxb, Dxf*Tepx^-1*Dyb, Dxf*Tepx^-1*Dzb; ...
        Dyf*Tepy^-1*Dxb, Dyf*Tepy^-1*Dyb, Dyf*Tepy^-1*Dzb; ...
        Dzf*Tepz^-1*Dxb, Dzf*Tepz^-1*Dyb, Dzf*Tepz^-1*Dzb];
    WAccelScal = speye(3*(Nx*Ny*Nz))*TmuSuper^-1;
    Ao = Ch*TmuSuper^-1*Ce - omega^2*TepsSuper;
    A = Ch*TmuSuper^-1*Ce + s*WAccelScal*GradDivMan - omega^2*TepsSuper;
    Ascal = Ch*Ce + s*GradDivMan - omega^2*TepsSuper*mu0;
    To = Ch*TmuSuper^-1*Ce + s*WAccelScal*GradDivHom -omega^2*TepsSuper;
    
%     figure; 
%     spy(A); pause; 

%% Construct Source Matrix
    J = [Mx; My; Mz];
    b = -1i*omega*J; bo = b;
    JCorrection = s * (1i/omega) * WAccelScal*TepsSuper^-1*GradDivHom*J;
    JCorrectionT = s* (1i/omega) * WAccelScal*TepsSuper^-1*GradDivHom*J;
    b = b+JCorrection;
    bT = b + JCorrectionT;
%% Construct Solutions
    solution = (A\b);
    solo = Ao\bo;
    solohom = To\bT;
    solLength = length(solution);
    Ex = solution(1:solLength/3);
    Ey = solution(solLength/3+1:solLength*(2/3));
    Ez = solution(solLength*(2/3)+1: solLength);
    Ez = reshape(full(Ez), Nx, Ny,Nz);
    figure;
    imagesc(abs(Ez(:,:,2)));
%% SOLVER CODE ENDS HERE
% 
% %% Visualize Eigenvalue Distribution: this gets clostly as matrix size scales
% if(Nx <= 10)
%    as1 = eig(full(A));
%    as2 = eig(full(Ao));
%    figure;
%    hist(abs(as1));
%    figure;
%    hist(abs(as2))
% end


%% visualize solutions
figure;
plot(abs(solution))
hold on;
plot(abs(solo))
legend('with wonwoeks accelerator', 'normal')
%% Do some Solution Visualizing
%[Exc, Eyc, Ezc] = FDFD3D_SliceVisualization(Ex, Ey, Ez, N, xrange, yrange);



