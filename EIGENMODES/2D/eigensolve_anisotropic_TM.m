

function [Hz_modes, Ex_modes, Ey_modes, eigenvals,A] = eigensolve_anisotropic_TM(L0, wvlen, xrange, ...
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
    exx = eps0*eps_tensor{1,1};
    exy = eps0*eps_tensor{1,2};
    eyx = eps0*eps_tensor{2,1};
    eyy = eps0*eps_tensor{2,2};
    N = size(exx);
    
    %% Set up the permittivity and permeability in the domain.
    % this should be x and y too for exx and eyy respectively 
    exx_avg = bwdmean_w(exx, 'x');  % average eps for eps_x
    eyy_avg = bwdmean_w(eyy, 'y');  % average eps for eps_y
    exy =  interpolate_array(interpolate_array(exy, 'y', 'f'), 'x','b');  % average eps for eps_x
    eyx = interpolate_array(interpolate_array(eyx, 'x', 'f'), 'y','b');  % average eps for eps_y

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
    Eyx = reshape(eyx,M,1);
   
    %% create interpolation functions;
    Rxf = interpolate(N, 'x', 'f');
    Rxb = interpolate(N, 'x', 'b');
    Ryf = interpolate(N, 'y', 'f');
    Ryb = interpolate(N, 'y', 'b');
    
    Texx = spdiags(Exx,0,M,M); % creates an MxM matrix, which is the correct size,
    Teyy = spdiags(Eyy,0,M,M);
    %have to be careful...only invert nonzero exy or eyx
    Texy = Rxf*Ryb*spdiags(Exy, 0,M,M);
    Teyx = Ryf*Rxb*spdiags(Eyx, 0,M,M); %have large regions of the diagonal that are 0
%     Texy = spfun(@(x) 1./x, Texy);
%     Teyx = spfun(@(x) 1./x, Teyx);
    Texx_avg = spdiags(reshape(exx_avg,M,1),0,M,M);
    Teyy_avg = spdiags(reshape(eyy_avg,M,1),0,M,M);    
    %% final cell construction
    Tep = [Texx, Texy; Teyx, Teyy];
    %Tep_cell = {Texx, Texy; Teyx, Teyy};
    %% create the derivative oeprators w/ PML

    N = [Nx, Ny];
    dL = [dx dy]; % Remember, everything must be in SI units beforehand

    Dxf = createDws('x', 'f', dL, N); 
    Dyf = createDws('y', 'f', dL, N);
    Dyb = createDws('y', 'b', dL, N); 
    Dxb = createDws('x', 'b', dL, N); 
    Dxf = Sxf^-1*Dxf; Dyf = Syf^-1*Dyf;
    Dyb = Syb^-1*Dyb; Dxb = Sxb^-1*Dxb; 
    
    Cb = vertcat(Dyb, -Dxb);
    Cf = horzcat(-Dyf, Dxf);
    %% PROBLEM CONSTRUCTION
    disp('get operator');
    % TEP \ cb is the most expensive operation for this whole thing...
    % it is okay if the anisotropy is a small fraction of the domain...
    A = Cf*(Tep\Cb);
    disp('start eigensolve');
    
    %% eigensolver
    omega_est = 2*pi*c0/(wvlen);
    %[U,V] = eigs(A, neigs, 'smallestabs');
    %find eigenmodes near desired frequency
    [U,V] = eigs(A, neigs, omega_est^2*mu0);

    eigenvals = diag(V); %eigenvals solved are omega^2*mu0
    eigenvals = (1/mu0)*sqrt(eigenvals);
    
    %% process the eigenmodes
    for i = 1:neigs
        hz = U(:,i);
        Hz = reshape(hz, Nx, Ny);
        a1 = Dyb*hz; a2 = -Dxb*hz;
        a = [a1;a2];
        e_fields = Tep\a;
        ex = e_fields(1:M);
        ey = e_fields(M+1:end);
%         inv_ex = Texx+Texy*Teyy^-1*Teyx;
%         
%         inv_ey = Teyy-Teyx*Texx^-1*Texy;
%         ey= spfun(@(x) 1./x, inv_ey)*(a2 + Teyx*Texx^-1*a1);
%         ex = spfun(@(x) 1./x, inv_ex)*(a1-Texy*Teyy^-1*a2);
        Ex = reshape(ex, Nx, Ny);
        Ey = reshape(ey, Nx,Ny);
        
        Ex_modes{i} = Ex;
        Ey_modes{i} = Ey;
        Hz_modes{i} = Hz;

    end


end