
function [Hz_modes, Ex_modes, Ey_modes, eigenvals] = ...
    eigensolve_TM_dispersive_Kx(L0, omega, xrange, ...
    yrange, eps_r, Npml, neigs, kx_guess)
   
    %EIGENSOLVE_TM Summary of this function goes here
    %   Detailed explanation goes here
    % we solve for KX, Ky is implicitly assumed to be 0.
    % neigs: number of eigenvalues
    % omega is the INPUT
    % Q(a) = Ma^2 + Da + K
    
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

    I = speye(M);
    Z = zeros(M,M,'like',I);    
    A = (1/mu0)*(Dxf*Tex^-1*Dxb+Dyf*Tey^-1*Dyb) + omega^2*I; % 1, with PML, this is not necessarily symmetric
%   C = mu0^-1*Te_unavged^-1;                                % lambda^2: DIAGONAL MATRIX
    C = -mu0^-1*Tex^-1;
    %B = -2*(mu0^-1)*1i*Vxf*Tex^-1*Dxb;                          % lambda    
    B = (mu0^-1)*(1i*(Dxf*Tex^-1+Tex^-1*Dxb));
                          % lambda    
    LHS = [Z , -C; 
           I ,  Z];
    RHS = [A, B; 
           Z, I ];
    NL = size(LHS);
    
%     A = -(Dxf*Tex^-1*Dxb - (KX*Te_unavged^-1*KX) + ...
%          2*1i*KX * Vxf*Tex^-1*Dxb +...
%          Dyf*Tey^-1*Dyb)/mu0; %

        
    %% eigensolver

    disp('start eigensolve');
    %[U,V] = eigs(A, neigs, 'smallestabs');
    %find eigenmodes near desired frequency
    tic
    [U,V] = eigs(RHS,LHS,neigs,kx_guess); %polyeigs is actually not suitable because it will solve all eigens
    toc
    %we need to linearize the eigenproblem ourselves

    

    eigenvals = diag(V); %eigenvals solved are omega^2*mu0
    Hz_modes = cell(1);
    Ex_modes = cell(1);
    Ey_modes = cell(1);
    uz_bloch_modes = cell(1);
    ux_bloch_modes = cell(1);
    uy_bloch_modes = cell(1);

    %% process the eigenmodes
    
    x = linspace(xrange(1), xrange(2), N(1))+diff(xrange)/2; 
    %the fact that xrange is positive and negative in x will be a problem
    % when complex wavevectors are involved.
    
    x = repmat(x.', [N(2), 1]); 
    for i = 1:neigs
        Kx = V(i,i);
        hz = U(1:round(NL/2),i);%.*exp(-1i*real(Kx)*x);%.*exp(-1i*Kx*x);
        Hz = reshape(hz, Nx,Ny);
        uz_bloch_modes{i} = reshape(U(1:round(NL/2),i), Nx,Ny);
        ex = (1/(1i*omega))*Tey\(Dyb*hz);
        ey = -(1/(1i*omega))*Tex\(Dxb*hz);
        Ex = reshape(ex, Nx, Ny);
        Ey = reshape(ey, Nx,Ny);
        Hz_modes{i} = Hz;
        Ex_modes{i} = Ex;
        Ey_modes{i} = Ey;

    end
    
    
end

