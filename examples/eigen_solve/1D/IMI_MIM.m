close all
clear
L0 = 1e-6;

Nx = 500; 
eps0 = 8.854e-12*L0;
eps_diel = 16;
eps_r = eps_diel*ones(Nx,1);
c0 = 3e8;

n = 13;
lattice_constant = 1; 
Lx = lattice_constant; % microns

dx = Lx/Nx; % in microns;
xrange = [-Lx,Lx]/2;
x = linspace(xrange(1), xrange(2), Nx);
yrange = [0,0];

%omega units must be units of L0

%% drude
omega_p = 0.72*pi*1e15;
gamma = 5.5e12;
ky_eigs_store = []; ky_eigs_2 = [];
epsilon_tracker = [];
omega_scan = linspace(1e13, 1.2*omega_p, 100);
c = 0;
Npml = 50;
for omega = linspace(1e13, 1.2*omega_p, 100)
    eps_drude = 1; %1-omega_p^2/(omega^2-1i*gamma*omega);
    epsilon_tracker = [epsilon_tracker, eps_drude];
    eps_r(.4*Nx:.6*Nx) =  eps_drude;
    wvlen = 2*pi*c0/omega*1e6;
    [Hz, ky_eigs_c,A, Dxf, Dxb] = FDFD_1D_Ky_eigensolve(L0, dx, eps_r, omega, Nx, n, Npml);
    if(c<5 && c> 2)
       figure();
       for i = 1:n
            plot(x, real(Hz{i}))
            hold on;
       end
       drawnow();
    end%% parameter setup


    ky_eigs_store = [ky_eigs_store,  ky_eigs_c];
    c=c+1;
    
end

figure();
for i = 1:2
    scatter(real(ky_eigs_store(i,:)), omega_scan./(2*pi*c0)*1e-6, '.b')
    hold on;
    scatter(abs(imag(ky_eigs_store(i,:))), omega_scan./(2*pi*c0)*1e-6, '.r')
end

ky_space = linspace(0,60, Nx);
hold on;
omega_med = c0*ky_space/sqrt(eps_diel)/((2*pi*c0));
omega_light = c0*ky_space/((2*pi*c0));

%% plot various light lines
plot(ky_space, omega_med)
plot(ky_space, omega_light);

%% get the plasmon cutoffs
line([0, max(ky_space)], [omega_p, omega_p]/(2*pi*c0)*1e-6)
line([0, max(ky_space)], [omega_p, omega_p]/(2*pi*c0)*1e-6/sqrt(eps_diel+1))


ylabel('frequency/(2*pi*c0)')
xlabel('ky (\mu m^-1)')
title('FDFD IMI even parity band')
legend('real', 'imaginary')


