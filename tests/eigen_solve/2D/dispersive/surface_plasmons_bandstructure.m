%% test script for TM Ex Ey eigensolve

close all
clear

% for surface plasmons, you would want a PEC (no propagating waves from
% spps)

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
mu0 = 4*pi*1e-7*L0;
eps0 = 8.854e-12*L0;
c0 = 3e8;
xrange = 0.2*[-1,1];  % x boundaries in L0
yrange = 1*[-1,1];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [100 150];  % [Nx Ny]
Npml = 1*[0 10];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);

wvlen = 4;
%% Set up the permittivity.
omega = 2*pi*c0/(wvlen*L0);
omega_p =  0.72*pi*1e15;
gamma = 5.5e12;
epsilon_metal = 1-omega_p^2/(omega^2-1i*gamma*omega);
epsilon = ones(N);
x = 1:N(1);
y = 1:N(2);
[xx, yy] = meshgrid(x,y);
half_ny = 14;
xbounds = [-0.1, 0.1];
ybounds = [-half_ny*dL(2), half_ny*dL(2)];
epsilon(:, 1:cy) = epsilon_metal;
figure();
visreal(epsilon, xrange, yrange);
drawnow();

%% eigensolve
neigs = 60;
kx_guess = 0.5*pi/L(1);
eps_r=epsilon;
%% =========================================================================
%% =========================================================================
    eps_0 = 8.85*10^-12*L0;
    mu_0 = 4*pi*10^-7*L0; 
    eps0 = eps_0;  % vacuum permittivity
    mu0 = mu_0;  % vacuum permeability in
    c0 = 1/sqrt(eps0*mu0);  % speed of light in 
    %% unravel the epsilon tensor;
    N = size(eps_r);
    
    %% Set up the permittivity and permeability in the domain.
    % the average should be x and y...but in what case do we not want that?
    exx = bwdmean_w(eps_0*eps_r, 'x');  % average eps for eps_x
    eyy = bwdmean_w(eps_0*eps_r, 'y');  % average eps for eps_y

    %% Set up number of cells
    xmin = xrange(1); xmax = xrange(2);
    ymin = yrange(1); ymax = yrange(2);
    Nx = N(1); dx = (xmax-xmin)/Nx;
    Ny = N(2); dy = (ymax-ymin)/Ny;
    
    M = prod([Nx, Ny]); %total number of cells
    L = [diff(xrange), diff(yrange)];
    %% Set up the Split coordinate PML
    Nx_pml = Npml(1); Ny_pml = Npml(2);
    lnR = -12; m = 3.5;
    
    sxf = create_sfactor(xrange,'f',omega,eps_0,mu_0,Nx,Nx_pml, lnR, m);
    syf = create_sfactor(yrange,'f', omega,eps_0,mu_0,Ny,Ny_pml, lnR, m);
    sxb = create_sfactor(xrange, 'b', omega,eps_0,mu_0, Nx, Nx_pml, lnR, m);
    syb = create_sfactor(yrange,'b', omega,eps_0,mu_0,Ny,Ny_pml, lnR, m);

    % now we create the matrix (i.e. repeat sxf Ny times repeat Syf Nx times)
    [Sxf, Syf] = ndgrid(sxf, syf);
    [Sxb, Syb] = ndgrid(sxb, syb);

    %Sxf(:) converts from n x n t0 n^2 x 1
    Sxf=spdiags(Sxf(:),0,M,M);
    Sxb=spdiags(Sxb(:),0,M,M);
    Syf=spdiags(Syf(:),0,M,M);
    Syb=spdiags(Syb(:),0,M,M);

    
    %% Create the dielectric and permeability arrays (ex, ey, muz)
    Tex = spdiags(reshape(exx, M,1),0,M,M);
    Tey = spdiags(reshape(eyy, M,1),0,M,M);
    Te_unavged = spdiags(eps0*eps_r(:), 0,M,M);
    %% create the derivative oeprators w/ PML
    N = [Nx, Ny];
    dL = [dx dy]; % Remember, everything must be in SI units beforehand
   
%     Dxf = createDws_bloch('x', 'f', dL, N, K_vec, L); 
%     Dyf = createDws_bloch('y', 'f', dL, N, K_vec, L);
%     Dyb = createDws_bloch('y', 'b', dL, N, K_vec, L); 
%     Dxb = createDws_bloch('x', 'b', dL, N, K_vec, L); 

    Dxf = createDws('x', 'f', dL, N); 
    Dyf = createDws('y', 'f', dL, N);
    Dyb = createDws('y', 'b', dL, N); 
    Dxb = createDws('x', 'b', dL, N); 
    Vxf = abs(createDws('x', 'f', [2 2], N)); 

    Dxf = Sxf\Dxf; Dyf = Syf\Dyf;
    Dyb = Syb\Dyb; Dxb = Sxb\Dxb; 
    
    %% PROBLEM CONSTRUCTION
    % assume H_zk(x) = e^{ikx}*(u_k(x)), and perform derivatives

    disp('get operator');
    I = speye(M);
    Z = zeros(M,M,'like',I);    
    %this appears to be a largely incorrect, though it's not easy to see
    %why. obviously, the dws_bloch does not explicitly know about the bloch
    % function form, but that's fine...we've directly solved that

    K = (1/mu0)*(Dxf*Tex^-1*Dxb+Dyf*Tey^-1*Dyb) + omega^2*I; % 1, with PML, this is not necessarily symmetric
    M = mu0^-1*Te_unavged^-1;                                % lambda^2: DIAGONAL MATRIX
    %D = -2*(mu0^-1)*1i*Vxf*Tex^-1*Dxb;% lambda    
    %apparently, these are not the same...
    D = -2*(mu0^-1)*Vxf*(1i*(Dxf*Tex^-1 + Tex^-1*Dxb));
    G = [M,Z; Z,I];
    C = [D,K; -I,Z];
    NL = size(G);
    
%     A = -(Dxf*Tex^-1*Dxb - (KX*Te_unavged^-1*KX) + ...
%          2*1i*KX * Vxf*Tex^-1*Dxb +...
%          Dyf*Tey^-1*Dyb)/mu0; %

        
    %% eigensolver

    disp('start eigensolve');
    %[U,V] = eigs(A, neigs, 'smallestabs');
    %find eigenmodes near desired frequency
    %Av = lambdaBv
    [U,V] = eigs(C,G,neigs,kx_guess); %polyeigs is actually not suitable because it will solve all eigens
    %we need to linearize the eigenproblem ourselves
    

    eigenvals = diag(V); %eigenvals solved are omega^2*mu0
    Hz_modes = cell(1);
    Ex_modes = cell(1);
    Ey_modes = cell(1);
    %% process the eigenmodes
    
    x = linspace(xrange(1), xrange(2), N(1)); 
    x = repmat(x.', [N(2), 1]); 
    for i = 1:neigs
        Kx = V(i,i);
        hz = U(1:round(NL/2),i); %.*exp(-1i*Kx*x);
        Hz = reshape(hz, Nx,Ny);
        ex = (1/(1i*omega))*Tey\(Dyb*hz);
        ey = -(1/(1i*omega))*Tex\(Dxb*hz);
        Ex = reshape(ex, Nx, Ny);
        Ey = reshape(ey, Nx,Ny);
        Hz_modes{i} = Hz;
        Ex_modes{i} = Ex;
        Ey_modes{i} = Ey;

    end
    

kx_eigs = eigenvals

%% =========================================================================
%% =========================================================================
[filtered_modes, filtered_k] = ...
<<<<<<< HEAD:tests/eigen_solve/2D/dispersive/surface_plasmons_bandstructure.m
    mode_filtering(Hz_modes, kx_eigs, epsilon, xbounds, ybounds, L, Npml);
=======
    mode_filtering(Hz_modes, kx_eigs, epsilon, xlim, ylim, L, Npml);
>>>>>>> dispersive_eigensolver:tests/eigen_solve/2D/dispersive/test_dispersive_2D_TM_eigensolve.m
filtered_modes = Hz_modes; filtered_k = kx_eigs;
for i = 1:length(filtered_k)
    figure();
    Kx = filtered_k(i);
    visreal(filtered_modes{i}, xrange, yrange);
    title(strcat(num2str(i), ', ', num2str(real(Kx)/(2*pi)*diff(xrange))));  
end
