%% test script for TM Ex Ey eigensolve

close all
clear

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
c0 = 3e8;
xrange = [-1,1];  % x boundaries in L0
yrange = [-1,1];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [200 200];  % [Nx Ny]
Npml = 1*[15 15];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);

%% Set up the permittivity.
wvlen_guess = 2;

epsilon = ones(N);
x = 1:N(1);
y = 1:N(2);
[xx, yy] = meshgrid(x,y);
epsilon((xx-cx).^2+(yy-cy).^2 < 50^2)=12;
epsilon((xx-cx).^2+(yy-cy).^2 < 30^2) = 1;
figure();
imagesc(epsilon);
drawnow();

%% eigensolve
neigs = 20;

[Ez, Hx, Hy, eigenvals] = eigensolve_TE(L0, wvlen_guess, xrange, ...
    yrange, epsilon, Npml, neigs);

for i = 1:neigs
   figure();
   subplot(121)
   visreal(Ez{i}, xrange, yrange);
   subplot(122)
   visreal(Hx{i}, xrange, yrange);
   title(eigenvals(i));

end
