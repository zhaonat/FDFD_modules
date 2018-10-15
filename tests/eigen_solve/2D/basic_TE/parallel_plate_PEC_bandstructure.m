%% test script for TM Ex Ey eigensolve

close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;
xrange = 0.05*[-1,1];  % x boundaries in L0
yrange = 0.5*[-1,1];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [80 120];  % [Nx Ny]
Npml = 0*[10 10];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);

%% Set up the permittivity.

epsilon = ones(N);
x = 1:N(1);
y = 1:N(2);

figure();
imagesc(epsilon);
drawnow();

neigs = 10;
kx_scan = linspace(-pi/L(1),pi/L(1), 30);
figure();
wvlen = 1;
for Kx = kx_scan;
    %% =======================================================================
    %% eigensolve
    %% =======================================================================
    K_vec = [Kx,0];

    %EIGENSOLVE_TM Summary of this function goes here
    %   Detailed explanation goes here
    eps_iso = epsilon;
    eps_0 = 8.85*10^-12*L0;
    mu_0 = 4*pi*10^-7*L0; 
    eps0 = eps_0;  % vacuum permittivity
    mu0 = mu_0;  % vacuum permeability in
    c0 = 1/sqrt(eps0*mu0);  % speed of light in 
    omega = 2*pi*c0/(wvlen);  % angular frequency in rad/sec

    %% unravel the epsilon tensor;
    N = size(eps_iso);

    %% Set up the permittivity and permeability in the domain.
    ezz = bwdmean_w(eps0*eps_iso, 'z');  % average eps for eps_y

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
    Tez = spdiags(reshape(ezz, M,1), 0, M,M);


    %% create the derivative oeprators w/ PML
    N = [Nx, Ny];
    dL = [dx dy]; % Remember, everything must be in SI units beforehand

    Dxf = createDws_bloch('x', 'f', dL, N, K_vec, L); 
    Dyf = createDws_bloch('y', 'f', dL, N, K_vec, L);
    Dyb = createDws_bloch('y', 'b', dL, N, K_vec, L); 
    Dxb = createDws_bloch('x', 'b', dL, N, K_vec, L); 
    % Dxf = createDws('x', 'f', dL, N); 
    % Dyf = createDws('y', 'f', dL, N);
    % Dyb = createDws('y', 'b', dL, N); 
    % Dxb = createDws('x', 'b', dL, N); 
    Dxf = Sxf^-1*Dxf; Dyf = Syf^-1*Dyf;
    Dyb = Syb^-1*Dyb; Dxb = Sxb^-1*Dxb; 

    %% PROBLEM CONSTRUCTION
    disp('get operator');

    %% create PEC mask
    M = prod(N);
    mask = ones(N);
    xn = 1:N(1);
    yn = 1:N(2);
    [Xn,Yn] = meshgrid(xn,yn);
    Xn = Xn.'; Yn = Yn.';
    mask(Yn == 1) = 0;
    mask(Yn == N(2)) =0;
    PEC_mask = spdiags(mask(:),0,M,M);

    A = -(1/mu0)*PEC_mask*(Dxb*Dxf + Dyb*Dyf)*PEC_mask; %
    B = Tez;

    disp('start eigensolve');
    %% eigensolver
    omega_est = 2*pi*c0/(wvlen);
    %[U,V] = eigs(A, B, neigs, 'smallestabs');
    %find eigenmodes near desired frequency
    [U,V] = eigs(A, B, neigs, omega_est^2);

    eigenvals = sqrt(diag(V)); %eigenvals solved are omega^2*mu0

    Ez_modes= cell(1);
    Hx_modes = cell(1);
    Hy_modes =cell(1);
    %% process the eigenmodes
    for i = 1:neigs
        ez = U(:,i);
        omega = V(i,i);
        Ez = reshape(ez, Nx,Ny);
        hx = (1/(1i*mu0*omega))*Dyf*ez;
        hy = -(1/(1i*mu0*omega))*Dxf*ez;
        Hx = reshape(hx, Nx, Ny);
        Hy = reshape(hy, Nx,Ny);
        Ez_modes{i} = Ez;
        Hx_modes{i} = Hx;
        Hy_modes{i} = Hy;

    end


    %% =======================================================================
    %% end eigensolve
    %% =======================================================================
    for i = 1:length(eigenvals)
        if(real(eigenvals(i)) == 0)
           continue; 
        end
        plot(Kx*L(1)/pi, real(eigenvals(i)), '.')
        hold on;
        drawnow();
    end
end

for i = 1:neigs
%    if(sum(isnan(Hx_modes{i}))>0)
%         continue;
%    end
    if(real(eigenvals(i)) <1e10)
       continue; 
    end
   figure();
   subplot(121)
   visreal(Ez_modes{i}, xrange, yrange);
   subplot(122)
   visreal(Hy_modes{i}, xrange, yrange);
   title(eigenvals(i));

end