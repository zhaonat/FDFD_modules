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
N = [80 120];  % [Nx Ny]
Npml = 1*[0,20];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
epsilon = ones(N);
thickness = 1;
d= thickness;
within_fabry = @(x,y) y < d/2 & y>-d/2;
epsilon = assign_val(epsilon, xrange, yrange, within_fabry, 12);

%% PROBES
probe_ind_y_ref = N(2)-Npml(2)-10;

%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = [1, Npml(2)+10];  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(:, ind_src(2)) = 1;

c = 1
wvlen_scan = linspace(1,2,100);
for wvlen = wvlen_scan

    k0 = 2*pi/wvlen;
    omega = 2*pi*c0/wvlen*1e6;
    epsilon_diel = 12;
    Hz_fields = cell(1);
    S_cell = cell(1);

    %% Solve TE (photonic), TM (practical), depends on convention equations.
    tic
    [Hz{c}, Ex, Ey] = ...
        solveTE(L0, wvlen, xrange, yrange, epsilon, Mz, Npml);
    toc

    Hz_fields{c} = Hz;

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
        solveTE(L0, wvlen, xrange, yrange, ones(N), Mz, Npml);
    toc
    [Sxr, Syr] = poynting(Hzr{c}, Exr, Eyr);

    
    %% Reflection

    R_vec(c) = abs(sum(Sy(:, probe_ind_y_ref))) / abs(sum(Syr(:, probe_ind_y_ref)))

    c = c+1;
    %drawnow();
end

f = figure()
plot(wvlen_scan, R_vec)
hold on;
plot(wvlen_scan, 1-R_vec);
legend( 'transmission', 'reflection')




