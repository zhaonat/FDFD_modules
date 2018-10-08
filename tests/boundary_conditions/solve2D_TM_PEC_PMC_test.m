%% A PMC is just the same as the PEC
% but a PEC for the TM polarization is not trivial
% we are solving for the Hz field in the scalar equation, but the Ex field
% and Ey fields are the tangential fields and they are only tangential to
% certain boundaries (not all of them, so doing a mask seems wrong)

close all
clear all 

%% Parameter Setup
L0 = 1e-6;  % length unit: microns
wvlen = 1;  % wavelength in L0
xrange = 1.5*[-1 1];  % x boundaries in L0
yrange = 1.5*[-1 1];  % y boundaries in L0

%% ASSYMETRIC FIELD PATTERNS CAN RESULT IF THE SOURCE IS NOT PERFECTLY CENTERED
% (which IT WON'T BE EVEN IF YOU PUT IT AT THE CENTRAL NODE POINT
% it seems that the node of the center in the scalar Hz field formulation
% will never be perfectly in the center of the grid?? which is weird...
% espeically consideri

Nx = 600; Ny = Nx;
N = [Nx Ny];  % [Nx Ny]
Npml = 0*[10 10];  % [Nx_pml Ny_pml]
mu0 = 4*pi*10^-7*L0; 
eps0 = 8.85*10^-12*L0;
c0 = 1/sqrt(mu0*eps0);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
eps_r = ones(N);

%% Fancy Structure

% eps_r(Nx/2:Nx/2+2,:) = 12;
% for j = 1:Ny
%     if(mod(j,10) == 0 || mod(j,11) == 0 || mod(j, 12) == 0)
%        eps_r(Nx/2:Nx/2+10,j) = 0; 
%     end
% end

%% Set up the 1agnetic current source density.
Mz = zeros(N);
ind_src = [N(1)/2 N(2)/2];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;
%a line source is also useful
%Jz(ind_src(1), :)) = 1;

%% SOLVER CODE STARTS HERE
    %normal SI parameters
    c0 = 1/sqrt(eps0*mu0);  % speed of light in 
    N = size(eps_r);  % [Nx Ny] THIS IS THE POINT WHERE THE GRID SIZE IS DETERMINED
    omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec

    %% Set up the permittivity and permeability in the domain.
    % bwdmean does nearest neighbor averaging (smoothes out stuff)

    eps_z = bwdmean_w(eps0 *eps_r, 'z'); %doesn't do anything in 2d
    eps_x = bwdmean_w(eps0 *eps_r, 'x');
    eps_y = bwdmean_w(eps0 *eps_r, 'y');

    %these are fully dense matrices...


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
    sxf = create_sfactor_mine(xrange,'f',omega,eps0,mu0,Nwx,Nx_pml);
    syf = create_sfactor_mine(yrange,'f', omega,eps0,mu0,Nwy,Ny_pml);
    sxb = create_sfactor_mine(xrange, 'b', omega,eps0,mu0, Nwx, Nx_pml);
    syb = create_sfactor_mine(yrange,'b', omega,eps0,mu0,Nwy,Ny_pml);

    % now we create the matrix (i.e. repeat sxf Ny times repeat Syf Nx times)
    [Sxf, Syf] = ndgrid(sxf, syf);
    [Sxb, Syb] = ndgrid(sxb, syb);

    %Sxf(:) converts from n x n t0 n^2 x 1
    Sxf=spdiags(Sxf(:),0,M,M);
    Sxb=spdiags(Sxb(:),0,M,M);
    Syf=spdiags(Syf(:),0,M,M);
    Syb=spdiags(Syb(:),0,M,M);


    %% Create the dielectric and permeability arrays (ex, ey, muz)
    %create a diagonal block matrix of ep and mu...
    epzList = reshape(eps_z,M,1); 
    epxList = reshape(eps_z,M,1); 
    epyList = reshape(eps_z,M,1); 

    Tepz = spdiags(epzList,0,M,M); % creates an MxM matrix, which is the correct size,
    Tepx = spdiags(epxList,0,M,M);
    Tepy= spdiags(epyList,0,M,M);
    %the M entries in epsList is put on the diagonals
    Tmz = mu0*speye(M); %in most cases, permeability is that of free-space
    Tmy = Tmz; Tmx = Tmz;

    %% Create Magnetic vector Mz (source profile determined by Mz input)
    % dimension = M*1
    Mz = reshape(Mz,M,1);
    Mz = sparse(Mz);

    %% create the derivative oeprators w/ PML

    N = [Nx, Ny];
    dL = [dx dy]; % Remember, everything must be in SI units beforehand

    Dxf = createDws('x', 'f', dL, N); 
    Dyf = createDws('y', 'f', dL, N);
    Dyb = createDws('y', 'b', dL, N); 
    Dxb = createDws_dirichlet_2D('x', 'b', dL, N); 
    Dxf_pml = Sxf^-1*Dxf; 
    Dyf_pml = Syf^-1*Dyf;
    Dyb_pml = Syb^-1*Dyb; 
    Dxb_pml = Sxb^-1*Dxb;
    
    %% construct PEC mask
    xn = 1:N(1);
    yn = 1:N(2);
    [Xn,Yn] = meshgrid(xn,yn);
    Xn = Xn.'; Yn = Yn.';
    maskx = ones(N);
    maskx(Xn == 1) = 0;
    maskx(Xn == N(1)) =0;
    
    %% ==================================================================
    %% TEST THE MASK, PEC or PMC
    %% ==================================================================

    % right now, if we wrap the entire grid in a PEC, the field pattern is
    % not symmetric...for sufficiently large domain size to wavelength
    % consequence of dispersion?
    masky = ones(N);
    masky(Yn == 1) = 0;
    masky(Yn == N(2)) =0;
    
    PMC_mask_x = spdiags(maskx(:), 0, M,M);
    PMC_mask_y = spdiags(masky(:), 0, M,M);
    
    % now I want a PEC mask for this polarization... perhaps we need to
    % reference the original FDFD equations...
    
    %% ==================================================================
    %% ==================================================================


    %% Construct the matrix A, everything is in 2D
    % this is the TE mode...
    A =  PMC_mask_y*PMC_mask_x*(Dxb_pml*(Tepx^-1)*Dxf_pml+ ...
        Dyb_pml*(Tepy^-1)*Dyf_pml)*PMC_mask_x*PMC_mask_y+ omega^2*Tmz;
    %A = PEC_mask_y*PEC_mask_x*A*PEC_mask_x*PEC_mask_y; 
    A0 = (Dxb_pml*(Tepx^-1)*Dxf_pml + Dyb_pml*(Tepy^-1)*Dyf_pml) + omega^2*Tmz;

    %% construct the matrix b, everything is in 2D
    b = 1i*omega*Mz;

    %% solve system
     tic
     % %% Solve the equation.
     if all(b==0)
        ez = zeros(size(b));
     else
       %hz = A\b;
        ez = A\b;
     end
     toc
     Ez = reshape(ez, N);

     %% now solve for Ex and Ey
     hx = -1/(1i*omega)*(Tmx^-1*Dyf)*ez;
     hy = (Tmy^-1*Dxf)*ez*(1/(1i*omega));
     Hy = reshape(hy,N);
     Hx = reshape(hx,N);

%% SOLVER CODE ENDS

 
 figure;
 visabs((Ez), xrange, yrange)
 figure;
 visabs((Hy), xrange, yrange)
 figure;
 visabs((Hx), xrange, yrange)
 
 moviereal(Ez, xrange, yrange)
 
 