close all
clear

%% general tests to verify an eigensolver is working (aside from looking at the
% field patterns

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
mu0 = 4*pi*1e-7*L0;
eps0 = 8.854e-12*L0;
c0 = 3e8;
xrange = 0.1*[-1,1];  % x boundaries in L0
yrange = 0.5*[-1,1];  % y boundaries in L0
L = [diff(xrange), diff(yrange)];
N = [100 250];  % [Nx Ny]
Npml = 1*[0 10];  % [Nx_pml Ny_pml]

[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML
Nx = N(1); Ny = N(2);
cx = round(Nx/2); cy = round(Ny/2);
wvlen = 2.8;
%% Set up the permittivity.
omega = 2*pi*c0/(wvlen*L0);
omega_p =  0.72*pi*1e15;
gamma = 5.5e12;
epsilon_metal = 1-omega_p^2/(omega^2-1i*gamma*omega);
epsilon = ones(N);
x = linspace(xrange(1), xrange(2), Nx);
y = linspace(yrange(1), yrange(2), Ny);
[xx, yy] = meshgrid(x,y); 
xx = xx.'; yy = yy.';
half_ny = 20;
xlim = [-0.1, 0.1];
ylim = [-half_ny*dL(2), half_ny*dL(2)];
epsilon(:,cy-half_ny:cy+half_ny) = 16;
epsilon(cx-10:cx+10, cy-half_ny:cy+half_ny) = epsilon_metal;
figure();
visreal(epsilon, xrange, yrange);
drawnow();

%% eigensolve
neigs = 30;
kx_guess = 0*pi/L(1);
[Hz_modes, Ex_modes, Ey_modes, kx_eigs, LHS, RHS] = ...
    eigensolve_TM_dispersive_Kx(L0, omega, xrange, yrange, epsilon, Npml, neigs, kx_guess);

filtered_k = kx_eigs;
filtered_modes = Hz_modes;
for i = 1:length(filtered_k)
    figure();
    Kx = filtered_k(i);%small Kx, but these modes should not be asymmetric
    eigenVec = filtered_modes{i}(:);
    
    expandedVec = [eigenVec; Kx*eigenVec];
    plot(real(RHS*expandedVec));
    hold on;
    plot(real(Kx*LHS*expandedVec), 'o');
    
end
drawnow();

%% compare to the non-dispersive solver
% neigs = 60;
% [Ez, Hx, Hy, omega_eigs] = solveTM_BlochX(L0, wvlen, xrange, yrange, epsilon, kx_guess, Npml, neigs);
% [filtered_modes_check, filtered_omega] = ...
%     mode_filtering(Ez, omega_eigs, epsilon, xlim, ylim, L, Npml);
% 
% for i = 1:length(omega_eigs)
%     if(abs(real(omega_eigs(i)))/omega > 0.1 && abs(real(omega_eigs(i)))/omega <6 )
%         figure();
% 
%         visreal(Ez{i}, xrange, yrange);
%         title(strcat(num2str(i), ', ', num2str(omega_eigs(i))));
%     end
% end