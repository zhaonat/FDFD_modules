%% example shows how to calculate reflection and transmission
% using the typical example of a fabry perot

close all
clear
%% potential applications
% cavity mirrors: multi-frequency, lasers
% oled mirrors

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;
xrange = 0.4*[-1 1];  % x boundaries in L0
yrange = [-2 2];  % y boundaries in L0
L = [diff(xrange), diff(yrange)]
N = [80 220];  % [Nx Ny]
Npml = 1*[0,20];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
x_vec = linspace(xrange(1), xrange(2),N(1));
Nx = N(1); Ny = N(2);
epsilon = ones(N);
thickness = 1;
d= thickness;
within_fabry = @(x,y) y < d/2 & y>-d/2;
epsilon = assign_val(epsilon, xrange, yrange, within_fabry, 12);

%% PROBES
probe_ind_y = Npml(2)+5;
probe_ind_y_ref = N(2)-Npml(2)-10;

%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = [1, Npml(2)+30];  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(:, ind_src(2)) = 1;

figure(); imagesc(abs(Mz.'));
hold on;
line([1, N(1)],[probe_ind_y, probe_ind_y])
line([1, N(1)],[probe_ind_y_ref, probe_ind_y_ref])


c = 1
degrees = pi/180;
wvlen_scan = linspace(1,1.5,50);
theta = 8*degrees;

for wvlen = wvlen_scan

    k0 = 2*pi/wvlen;
    omega = 2*pi*c0/wvlen*1e6;
    epsilon_diel = 12;
    Hz_fields = cell(1);
    S_cell = cell(1);
    kx = 2*pi/wvlen*sin(theta);
    ky = 2*pi/wvlen*cos(theta);
    Kvec = [kx,ky];
    
    %% modify source
    Mz = zeros(N);
    Mz(:, ind_src(2)) = exp(-1i*(kx*x_vec.'));
    Mz(:, ind_src(2)-1) = Mz(:, ind_src(2)).*exp(-1i*(ky * dL(2)) -1i*pi);
    
    
    %% Solve TE (photonic), TM (practical), depends on convention equations.
    tic
    [Hz{c}, Ex, Ey] = ...
        solveTE_bloch(L0, wvlen, xrange, yrange, epsilon, Mz, Npml, Kvec);
    toc


    %% Visualize the solution.
    % figure;
    % visreal(Hz, xrange, yrange)

    %%
    [Sx,Sy] = poynting(Hz{c}, Ex, Ey);
    S_cell{c,1} =Sx;
    S_cell{c,2} = Sy;

    %%reference simulation
    tic
    [Hzr{c}, Exr, Eyr] = ...
        solveTE_bloch(L0, wvlen, xrange, yrange, ones(N), Mz, Npml,Kvec);
    toc
    [Sxr, Syr] = poynting(Hzr{c}, Exr, Eyr);

    
    %% Reflection

    R_vec(c) = abs(sum(Sy(:, probe_ind_y))) / abs(sum(Syr(:, probe_ind_y_ref)))

    c = c+1;
    %drawnow();
end

figure()
plot(wvlen_scan, R_vec)
hold on;
plot(wvlen_scan, 1-R_vec);
legend('R', 'T');



