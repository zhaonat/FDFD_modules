%% test script for TM Ex Ey eigensolve
%% getting plasmonic modes is really hard
close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;
xrange = [-1,1];  % x boundaries in L0
yrange = [-2,2];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [200 200];  % [Nx Ny]
Npml = 1*[0 20];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);

%% Set up the permittivity.
wvlen_guess = 0.5;

epsilon = ones(N);
x = 1:N(1);
y = 1:N(2);
[xx, yy] = meshgrid(x,y);
epsilon(:, 1:floor(Ny/2)) = -2-0.01i;
figure();
visreal(epsilon, xrange, yrange);
drawnow();

%% eigensolve
neigs = 20;
K_vec = [0,0];

[Hz, Ex, Ey, eigenvals] = eigensolve_TM(L0, wvlen_guess, xrange, ...
    yrange, epsilon, Npml, neigs, K_vec);

% for i = 1:neigs
%    figure();
%    subplot(121)
%    visreal(Hz{i}, xrange, yrange);
%    subplot(122)
%    visreal(Ex{i}, xrange, yrange);
%    title(eigenvals(i)/(2*pi*c0)*1e-6);
% 
% end

Kx = 0;
[Hzc, Exc, Eyc, eigenvalsc] = solveTE_BlochX(L0,...
    wvlen_guess, xrange, yrange, epsilon, Kx, Npml, neigs);

for i = 1:neigs
   figure();
   subplot(121)
   visreal(Hz{i}, xrange, yrange);
   subplot(122)
   visreal(Hzc{i}, xrange, yrange);
   title(eigenvals(i)/(2*pi*c0)*1e-6);

end
