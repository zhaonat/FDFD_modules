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

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
L = [diff(xrange), diff(yrange)];

%% construct nonuniform domain parameters

dx_scale = ones(1,N(1)); %start with some dx and scale each one pointwise
dy_scale = ones(1,N(2));

grading_N = 20; 
grade_up = logarithmic_grading(0.2, 1,grading_N);
grade_down = logarithmic_grading(1, 0.2,grading_N);
disp = 50;

dy_scale(disp:disp+length(grade_up)-1) = grade_down;
dy_scale(end-length(grade_up)-disp:end-disp-1) = grade_up;
dy_scale(disp+length(grade_up)-1:end-length(grade_up)-disp) = grade_down(end);
dx_scale(disp:disp+length(grade_up)-1) = grade_down;
dx_scale(end-length(grade_up)-disp:end-disp-1) = grade_up;
dx_scale(disp+length(grade_up)-1:end-length(grade_up)-disp) = grade_down(end);

dxvector = dx_scale*dL(1);
dyvector = dy_scale*dL(2);

%view the final grading 
figure();
subplot(121)
plot(dy_scale*dL(2)); subplot(122); plot(dx_scale*dL(1));
drawnow();

%convert xs to an xrange_array
xrange_array = cumsum(dx_scale)*dL(1);
yrange_array = cumsum(dy_scale)*dL(2);

figure();
plot(xrange_array);
hold on;
plot(yrange_array);

%% once we specify the grading, we can recalculate the xrange
Lx = sum(dx_scale*dL(1));
Ly = sum(dy_scale*dL(2));

xrange2 = 0.5*[-Lx, Lx];
yrange2 = 0.5*[-Ly, Ly];

%% test operator creator
% [Fsx, Fsy, Fsx_conj, Fsy_conj] = non_uniform_scaling(dx_scale, dy_scale);
% figure(); 
% subplot(221);
% plot(diag(Fsx));
% subplot(222);
% plot(diag(Fsy));
% subplot(223);
% plot(diag(Fsx_conj));
% subplot(224);
% plot(diag(Fsy_conj))
% 
% % dxf operators
% Dyf = createDws('y', 'f', dL, N);%*Fsx; 
% Dyft = createDws_nm('y', 'f', dyvector.', N);
% Dxft = createDws_nm('x','f', dxvector.', N);
% Dxf = createDws('x','f', dL, N);
% figure();
% plot(diag(Fsy^-1*Dyf), '.-g', 'markersize', 3); hold on; plot(diag(Dyft), '.-'); plot(diag(Dyf));
% figure();
% plot(diag(Dxft), '.-', 'markersize', 15);
% hold on;
% plot(diag(Fsx^-1*Dxf), '-mx')
% plot(diag(Dxf));
% %xlim([0,2*N(1)]);

%% Set up the permittivity.
wvlen = 0.75;
eps = ones(N);
%eps(50:150, 50:150) = 6;

%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = [round(Nx/2),round(Ny/2)];  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;

%% Solve TE (photonic), TM (practical), depends on convention equations.
% this is a bit weird: the grid dimension is Nx, Ny, same as uniform
% but it should be larger
% dr_corase is the same as the uniform, but dr_fine is not

tic
[Hz, Ex, Ey, A_nu,b_nu, Dxf, Dyf, Dxb, Dyb] = solveTE_nu(L0, wvlen, xrange, yrange,...
    eps, Mz, Npml, dx_scale, dy_scale);
toc
% figure();
% moviereal(Hz, xrange2, yrange2);

%% uniform
% since we use xrange and yrange the same as nonuniform...this should
% be the same exact system...
tic
[Hz_u, Ex_u, Ey_u,  A] = solveTE(L0, wvlen, xrange2, yrange2, eps, Mz, Npml);
toc
dlux = diff(xrange2)/N(1);
dluy = diff(yrange2)/N(2);
%
% figure();
% moviereal(Hz_u, xrange2, yrange2);

figure();
show_grid = 1;
visreal_nu(Hz, xrange_array, yrange_array, show_grid);
title('non-uniform')

figure();
visreal(Hz_u, xrange2, yrange2);
title('uniform')
figure();
visreal_nu(Hz, xrange_array, yrange_array, 0);
