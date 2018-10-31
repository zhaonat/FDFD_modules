close all
clear

%% testing the non-uniform mesh in real FDFD simulations


%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;

%% non-uniform grid specification
% 1) partition into coarse vs fine (vs pml)
% along a line, it looks like Ncoarse, Nfine, Ncoarse
% in this construction scheme, we actually do not set the full domain size
% in the initial case

%% a priori, we know what the physical parameter space looks like

Nfine = [40,100]; %specify nx, ny for each region
Ncoarse = [40,50];
Ntran =   [30,50];

%% LARGE DISCREPANCY IN DX AND DY CAN CUASE STRANGE EFFECTS!!!

% 2) specify the dx and dy of each region
dx1 = 0.002; dy1 = 0.02;
dx2 = 0.002; dy2 = 0.0025;
dfine = [dx2, dy2];
dcoarse = [dx1, dy1];
dtran = [0 ,0]; %nothing can go here since the transition has varying dims

% 3) stack the vectors
% drt does not have a value...
Nft = vertcat(Ncoarse, Ntran, Nfine, Ntran, Ncoarse);
drt = vertcat(dcoarse, dtran, dfine, dtran, dcoarse);

% scale is arbitrary, just take dcoarse;
dr_reference  = dcoarse;

% 4) construct scaling vectors from this information
[dx_scale, dy_scale] = generate_nonuniform_scaling(Nft, drt./dr_reference);

%% calculate Ntot and Ltot
N = sum(Nft);
Lx = sum(dr_reference(1)*dx_scale);
Ly = sum(dr_reference(2)*dy_scale);
L = [Lx, Ly];
xrange = 0.5*[-Lx, Lx];
yrange = 0.5*[-Ly, Ly];
xrange_array = cumsum(dr_reference(1)*dx_scale)-Lx/2;
yrange_array = cumsum(dr_reference(2)*dy_scale)-Ly/2;
Nx = N(1); Ny = N(2);
%% output is a dxscale...dyscale


%% PML specification
Npml = [0,20];

%% Set up the permittivity.
wvlen = 1.5;
omega = 2*pi*c0/wvlen*1e6;
eps = ones(N);
y_thickness = 0.3;

% x_center and y_center should be values ACTUALLY in xrange array
% it should be the users responsibility to do that...
x_center = 0;
y_center = -0.0;
fill_factor = 0.2;
epsilon_edge = 12;
omega_p = 0.72*pi*1e15;
gamma = 5.5e12;
epsilon_center= 1-omega_p^2/(omega^2-1i*omega*gamma);
N = size(eps);
Nx = N(1); Ny = N(2);
Lx = xrange_array(end)-xrange_array(1); 
Ly = yrange_array(end)-yrange_array(1);



%% multiple unit cells now
Nx = length(xrange_array); Ny = length(yrange_array);
eps_multi = ones(Nx, Ny);
num_cells = 2;

Lx = xrange_array(end)- xrange_array(1);
lattice_constant = Lx/num_cells;
%% now we have to stripe it, which we will do per unit cell
% this means, we need to locate AT LEAST, the COM of each unit cell
% which seems to be a bad idea...
for i = 0:num_cells-1
    %determine x_center of each unit cell
    x_center = xrange_array(1)+lattice_constant*i+lattice_constant/2;
    eps_multi = add_grating_nu(eps_multi, xrange_array, yrange_array, epsilon_edge, epsilon_center, ...
fill_factor, y_thickness, y_center, x_center, lattice_constant);

figure(); 
visreal_nu(eps_multi, xrange_array, yrange_array,0);

end
alpha =1;
[X,Y] = meshgrid(xrange_array,yrange_array);
figure();
imagesc(Y.' > 4*alpha*X.'.^2); colorbar;
%transpose...
eps_multi(Y.' > 4*alpha*X.'.^2 -0.1) = 1;
figure(); 
subplot(121)
visreal_nu(eps_multi, xrange_array, yrange_array,1);
subplot(122);
visreal_nu(eps_multi, xrange_array, yrange_array,0);

%% test the solveTE_nu
%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = [round(Nx/2),N(2) - Npml(2)-10];  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(:, ind_src(2)) = 1;


%% Solve TE (photonic), TM (practical), depends on convention equations.
tic
[Hz, Ex, Ey, A_nu,b_nu, Dxf_nu, Dyf_nu] = solveTE_nu(L0, wvlen, xrange, yrange,...
    eps_multi, Mz, Npml, dx_scale, dy_scale, dr_reference);
toc

figure();
visreal_nu(Hz, xrange_array, yrange_array, 0);
drawnow();



