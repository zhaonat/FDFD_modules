close all
clear

%% testing the non-uniform mesh in real FDFD simulations


%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;

xrange = [-2,2];  % x boundaries in L0
yrange = [-2.5,2.5];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = 2*[100,100];  % [Nx Ny]
Npml = 1*[20 20];  % [Nx_pml Ny_pml]
Nx = N(1); Ny = N(2);
epsilon_diel = 16;

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
L = [diff(xrange), diff(yrange)];

%% construct nonuniform domain parameters
xs = ones(1,N(1));
ys = ones(1,N(2));

grading_N = 30; alpha = 0.98;
grade_up = logarithmic_grading(grading_N, alpha, 'u');
grade_down = logarithmic_grading(grading_N, alpha, 'd');
disp = 30;
ys(disp:disp+length(grade_up)-1) = grade_up;
ys(end-length(grade_up)-disp:end-disp-1) = grade_down;
ys(disp+length(grade_up)-1:end-length(grade_up)-disp) = grade_up(end);

%view the final grading 
figure();
plot(ys);
drawnow();

%convert xs to an xrange_array
xrange_array = cumsum(xs)*dL(1);
yrange_array = cumsum(ys)*dL(2);

[Xs, Ys] = meshgrid(xrange_array, yrange_array);
%% once we specify the grading, we can recalculate the xrange
Lx = sum(xs*dL(1));
Ly = sum(ys*dL(2));
xrange2 = 0.5*[-Lx, Lx];
yrange2 = 0.5*[-Ly, Ly];

%% test operator creator
[Fsx, Fsy, Fsx_conj, Fsy_conj] = non_uniform_scaling(xs, ys);

%% Set up the permittivity.
wvlen = 1.5;
omega_p = 0.72*pi*1e15;
gamma = 5.5e12;
omega = 2*pi*c0/wvlen*1e6;
epsilon_metal =  1 - omega_p^2./(omega^2-1i*gamma*omega);
epsilon_diel = 16;

eps = ones(N);

%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = [round(Nx/2),round(Ny/2)];  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;

%% Solve TE (photonic), TM (practical), depends on convention equations.
tic
[Hz, Ex, Ey, A_nu,b_nu, Dxf_nu, Dyf_nu] = solveTE_nu(L0, wvlen, xrange2, yrange2, eps, Mz, Npml, xs, ys);
toc

% figure();
% moviereal(Hz, xrange2, yrange2);

%% no non-uniform
tic
[Hz_u, Ex_u, Ey_u,  A, A_mode, b, T_eps_x, T_eps_y, Dyb, Dxb] = solveTE(L0, wvlen, xrange2, yrange2, eps, Mz, Npml);
toc
%
% figure();
% moviereal(Hz_u, xrange2, yrange2);

figure();
visreal_nu(Hz, xrange_array, yrange_array);
figure();
visreal(Hz_u, xrange2, yrange2);