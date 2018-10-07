%% test script for TM Ex Ey eigensolve

% close all
% clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
mu0 = 4*pi*1e-7*L0;
eps0 = 8.854e-12*L0;
c0 = 1/sqrt(mu0*eps0);
num_cells = 2;
xrange = 0.1*num_cells*[-1,1];  % x boundaries in L0
yrange = [-2,2];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [100 150];  % [Nx Ny]
Npml = 1*[0 10];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);
wvlen = 2;
%% Set up the permittivity.
omega = 2*pi*c0/(wvlen);
omega_p =  0.72*pi*1e15;
gamma = 5.e12;
epsilon_metal = 1-omega_p^2/(omega^2-1i*gamma*omega);
epsilon_diel = 16; 
%epsilon_diel = epsilon_metal;
fill_factor = 0.2;
thickness = 0.4;
eps_r = ones(N);
y_grid_center = L(2)/2;
y_center = y_grid_center-0.8;
eps_r = hybrid_grating_multi_unit_cell_add(eps_r,num_cells, N, L, epsilon_diel,...
    epsilon_metal, fill_factor, thickness, y_center);
y_center = y_grid_center+0.8;
eps_r = hybrid_grating_multi_unit_cell_add(eps_r,num_cells, N, L, epsilon_diel,...
    epsilon_metal, fill_factor, thickness, y_center);
xlim = xrange;
ylim = [-1.5, 1.5];
figure();
visreal(eps_r, xrange, yrange);
drawnow()
%% eigensolve
neigs = 100;
kx_guess = 0*pi/L(1);
%% ========================================================================
%% ========================================================================
%% ========================================================================
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

%% ========================================================================
%% ========================================================================
%% ========================================================================
kx_eigs = diag(V);
pml_threshold = 1e-3;
[filtered_modes, filtered_k, mask] = ...
    mode_filtering(Hz_modes, kx_eigs, eps_r, xlim, ylim, L, Npml);
%filtered_k = kx_eigs; filtered_modes = Hz_modes;
for i = 1:length(filtered_k)
    figure();
    Kx = filtered_k(i);
    visreal(filtered_modes{i}, xrange, yrange);
    title(strcat(num2str(i), ', ', num2str(imag(Kx)/(2*pi)*diff(xrange))));  
end
% 
figure();
moviereal(filtered_modes{2}, xrange, yrange)

%% compare with bloch solver
% neigs = 50;
% [Ez, Hx, Hy, omega_eigs] = solveTM_BlochX(L0, wvlen, xrange, yrange, eps_r, kx_guess, Npml, neigs);
% [filtered_modes, filtered_omega, mask] = ...
%     mode_filtering(Ez, omega_eigs, eps_r, xlim, ylim, L, Npml);
% 
% for i = 1:length(filtered_omega)
%     if(abs(real(filtered_omega(i)))/omega > 0.5 )
%         figure();
% 
%         visreal(filtered_modes{i}, xrange, yrange);
%         title(filtered_omega(i));
%     end
% end
