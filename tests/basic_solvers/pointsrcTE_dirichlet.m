%ptsrc_TE with a dirichlet boundary condition encoded into DWS
clear all; close all; clc;

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
wvlen = 1.0*L0;  % wavelength in L0
xrange = L0*[-5 5];  % x boundaries in L0
yrange = L0*[-5 5];  % y boundaries in L0
N = [30 30];  % [Nx Ny]
Npml = [5 5];  % [Nx_pml Ny_pml] need to deal with the special case where both are 0

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
eps_air = ones(N);
%{
for i = 1:N
   for j = 1:N
       if(j >10 && j<30)
          eps_air(i,j)= 1.5;
       end
   end
end
%}

%% Set up the magnetic current source density.
Mz = zeros(N);
ind_src = [10 10];%ceil(N/2);  % (i,j) indices of the center cell; Nx, Ny should be odd
Mz(ind_src(1), ind_src(2)) = 1;

%% Solve TE equations.
[Hz, Ex, Ey, A, omega,b, Dxf, Dyf, Dxb, Dyb] = solveTE_dirichlet(L0, wvlen, xrange, yrange, eps_air, Mz, Npml);
%% Visualize the solution.
figure;
visabs(Hz, xrange, yrange)

%%
figure;
visabs(Ex, xrange, yrange - 0.5*dL(2));

figure;
visabs(Ey, xrange - 0.5*dL(1), yrange);

%% Show the movie of the oscillating field.
figure;
moviereal(Hz, xrange, yrange)

figure;
imagesc(log(abs(Hz)))



