close all
clear all
% simple test of the adjoint method for a single focal point
% things to be aware of
% what should we initialize epsilon as?
% what kind of algorithms should we use for doing gradient descent

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
eps0 = 8.854e-12*L0;
mu0 = 4*pi*1e-7*L0;
c0 = 1/sqrt(mu0*eps0);
num_cells = 2;
xrange = [-1 1];  % x boundaries in L0
yrange = [-1 1];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [100 100];  % [Nx Ny]
N0 = N; 
Npml = 1*[15 15];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
M = prod(N);

%% set up initial epsilon
eps_r= ones(N);


%% source placement
wvlen = 0.6;
omega = 2*pi*c0/wvlen;
source_coord = [0,-0.9];
[njx, njy] = coord_to_grid(source_coord, N, xrange, yrange);
Mz = zeros(N);
Mz(njx, njy) = 1;

%% define coordinate to optimize
opt_coord = [0, 0.9];
[nox, noy] = coord_to_grid(opt_coord, N, xrange, yrange);
eta = zeros(N);
eta(nox, noy) = 1; 
eta = eta(:);

%% define optimization mask
mask = zeros(N);
mask(40:end-40, 50:end-50) = 1;

%% ENTER Optimization loop
epochs = 100;

for t =1:epochs

    %% get source field
    tic
    [Ez, Hx, Hy, A, omega,b]  = ...
    solveTM(L0, wvlen, xrange, yrange, eps_r, Mz, Npml);
    toc
    u0 = Ez(:);

    %% objective
    objective_eval = conj(eta.'*u0)*(eta.'*u0);
    %% get adjoint field;
    grad_objective = -2*(eta.'*u0)*eta;
    u_adj = A\grad_objective;

    Ezu = reshape(u_adj, Nx,Ny);
    % derivative of A with respect to the parameter, which is eps_r(nox, noy);
    dA_deps = omega^2*speye(M);

    %% gradient for epsilon
    % problem... it appears I get the same value for everything...
    derivative_of_objective = full(-omega^2*(eps0/L0)*real(Ez.*Ezu));
    
    %% update epsilon
    alpha = 1e-3;
    alpha = alpha*0.99;
    new_eps = alpha*derivative_of_objective;
    
    %% enforce some reality conditions (eps_r < 12 for example)
    new_eps(new_eps> 12) = 12;
    eps_r = eps_r - mask.*new_eps;
    if(mod(t,4) == 0)
        figure(); 
        subplot(121)
        visabs(Ez, xrange, yrange);
        subplot(122);
        visreal(eps_r, xrange, yrange);
        caxis([-1,6])
        drawnow();
    end
    
end

