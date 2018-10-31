%% test script for TM Ex Ey eigensolve

close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
mu0 = 4*pi*1e-7*L0;
eps0 = 8.854e-12*L0;
c0 = 1/sqrt(mu0*eps0);
xrange = [-1,1];  % x boundaries in L0
yrange = [-2,2];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [100 150];  % [Nx Ny]
Npml = 1*[0 10];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);
wvlen = 1.5;
%% Set up the permittivity.
omega = 2*pi*c0/(wvlen);
omega_p =  0.72*pi*1e15;
gamma = 5.5e12;
epsilon_metal = 1 - omega_p^2/(omega^2-1i*gamma*omega);
epsilon_diel = 16; %epsilon_metal;
fill_factor = 0.2;
thickness = 0.4;
eps_r = ones(N);
num_cells = 1;
y_grid_center = L(2)/2;
y_center_1 = y_grid_center-0.8;
eps_r = hybrid_grating_multi_unit_cell_add(eps_r,num_cells, N, L, epsilon_diel,...
    epsilon_metal, fill_factor, thickness, y_center_1);
y_center_2 = y_grid_center+0.8;
eps_r = hybrid_grating_multi_unit_cell_add(eps_r,num_cells, N, L, epsilon_diel,...
    epsilon_metal, fill_factor, thickness, y_center_2);
xlim = xrange;
ylim = [-1.5, 1.5];
figure();
visreal(eps_r, xrange, yrange);
drawnow()
%% eigensolve
neigs = 50;
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
    lnR = -12; m = 2;
    
    sxf = create_sfactor(xrange, 'f', omega,eps_0,mu_0,Nx,Nx_pml, lnR, m);
    syf = create_sfactor(yrange, 'f', omega,eps_0,mu_0,Ny,Ny_pml, lnR, m);
    sxb = create_sfactor(xrange, 'b', omega,eps_0,mu_0, Nx, Nx_pml, lnR, m);
    syb = create_sfactor(yrange, 'b', omega,eps_0,mu_0,Ny,Ny_pml, lnR, m);

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
   

    Dxf = createDws('x', 'f', dL, N); 
    Dyf = createDws('y', 'f', dL, N);
    Dyb = createDws('y', 'b', dL, N); 
    Dxb = createDws('x', 'b', dL, N); 

    Dxf = Sxf\Dxf; Dyf = Syf\Dyf;
    Dyb = Syb\Dyb; Dxb = Sxb\Dxb; 
    
    %% PROBLEM CONSTRUCTION
    % assume H_zk(x) = e^{ikx}*(u_k(x)), and perform derivatives

    I = speye(M);
    Z = zeros(M,M,'like',I);    
    A = (mu0^-1)*Te_unavged^-1*(Dxb*Dxf + Dyb*Dyf) + omega^2*I; % 1, with PML, this is not necessarily symmetric
    C = -(mu0^-1)*Te_unavged^-1;                                % lambda^2: DIAGONAL MATRIX
    B = (mu0^-1)*Te_unavged^-1*1i*(Dxf+Dxb);                           % lambda    
    LHS = [Z , -C; 
           I ,  Z];
    RHS = [A, B; 
           Z, I ];
    NL = size(LHS);
    
        
    %% eigensolver

    %[U,V] = eigs(A, neigs, 'smallestabs');
    %find eigenmodes near desired frequency
    tic
    [U,V] = eigs(RHS,LHS,neigs,kx_guess); %polyeigs is actually not suitable because it will solve all eigens
    toc
    %we need to linearize the eigenproblem ourselves
    

    eigenvals = diag(V); %eigenvals solved are omega^2*mu0
    Ez_modes = cell(1);
    Hx_modes = cell(1);
    Hy_modes = cell(1);
    Ez_test= cell(1);
    blochz = cell(1);
    %% process the eigenmodes
    
    x = linspace(xrange(1), xrange(2), N(1))+diff(xrange)/2; 
    x = repmat(x.', [N(2), 1]); 
    for i = 1:neigs
        Kx = V(i,i);
        % the bloch part of the wavefunction is tricky if the medium is 
        % lossy...
        %in order to make things consistent, we have to separate realk,
        %imagk
        % also... fabry perot...modes that oscillate up and down...
        realk = real(Kx);
        imagk = imag(Kx);
        
        ez = U(1:round(NL/2),i).*exp(-1i*abs(real(Kx))*x); %note that the bloch part is tough
        ezk = U(round(NL/2)+1:end,i).*exp(1i*abs(realk)*x).*exp(-abs(imagk)*abs(x));
        Ezk = reshape(ezk, Nx,Ny);
        Ez = reshape(ez, Nx,Ny);
        hx = (1/(1i*omega))*Tey\(Dyb*ez);
        hy = -(1/(1i*omega))*Tex\(Dxb*ez);
        Hx = reshape(hx, Nx, Ny);
        Hy = reshape(hy, Nx,Ny);
        Ez_modes{i} = Ez;
        Hx_modes{i} = Hx;
        Hy_modes{i} = Hy;
        Ez_test{i} = Ezk; blochz{i} = reshape(U(1:round(NL/2),i), Nx,Ny);
        disp(c0^2*real(Kx)^2)
    end

%% ========================================================================
%% ========================================================================
%% ========================================================================
kx_eigs = diag(V);
[filtered_modes, filtered_k, mask] = ...
    mode_filtering(Ez_test, kx_eigs, eps_r, xlim, ylim, L, Npml);
filtered_k = kx_eigs; filtered_modes = Ez_test;

% a great service could be done if we could write an eigensolver which
% preferentially does not return the spurious modes (particularly from the
% PML)

% imaginary kx or real kx will fail to capture fabry perot type modes...
for i = 1:length(filtered_k)
    figure();
    Kx = filtered_k(i);
    subplot(131)
    visreal(filtered_modes{i}, xrange, yrange);
    title(strcat(num2str(i), ', ', num2str(real(Kx)/(2*pi)*diff(xrange))));  
    hold on;
%     line([xrange(1), xrange(2)],[y_center_1, y_center_1]);
%     line([xrange(1), xrange(2)],[y_center_2, y_center_2])
    subplot(132)
    visreal(blochz{i}, xrange, yrange);
    title(strcat(num2str(i), ', ', num2str(imag(Kx)/(2*pi)*diff(xrange)))); 
    subplot(133)
    visreal(Ez_modes{i}, xrange, yrange);
    hold on;

end
% 
figure();
moviereal(filtered_modes{2}, xrange, yrange)

%% compare with bloch solver
% neigs = 100;
% [Ez, Hx, Hy, omega_eigs] = solveTM_BlochX(L0, wvlen, xrange, yrange, eps_r, kx_guess, Npml, neigs);
% [filtered_modes_check, filtered_omega] = ...
%     mode_filtering(Ez, omega_eigs, eps_r, xlim, ylim, L, Npml);
% 
% for i = 1:length(omega_eigs)
%     if(abs(real(omega_eigs(i)))/omega > 0.6 && abs(real(omega_eigs(i)))/omega <1.4 )
%         figure();
% 
%         visreal(Ez{i}, xrange, yrange);
%         title(omega_eigs(i));
%     end
% end
