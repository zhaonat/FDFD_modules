close all
clear
%% potential applications
% cavity mirrors: multi-frequency, lasers
% oled mirrors

%% design principle: we want magnitude of dielectric equal to magnitude of metallic dielectric?

%% Set up the domain parameters.
L0 = 1e-6;  % length unit: microns
wvlen = 2.2;  % wavelength in L0 
c0 = 3e8;
xrange = 0.1*[-1 1];  % x boundaries in L0
yrange = 0.5*[-1 1];  % y boundaries in L0
a = diff(xrange);

thickness = diff(yrange);

L = [diff(xrange), diff(yrange)]
N = [150 210];  % [Nx Ny]
Npml = 1*[0 40];  % [Nx_pml Ny_pml]

Nx = N(1); Ny = N(2);
[xrange, yrange, N, dL, Lpml] = domain_with_pml(xrange, yrange, N, Npml);  % domain is expanded to include PML

%% Set up the permittivity.
omega_p = 0.72*pi*1e15;%3e15;
gamma = 5.5e12;
omega = 2*pi*c0/wvlen*1e6;
epsilon_diel = 16;
epsilon_metal =  1 - omega_p^2./(omega^2+1i*gamma*omega);

fill_factor = 0.2; %half metal, half dielectric

eps = hybrid_grating_grid(N, L, epsilon_diel, epsilon_metal, fill_factor, diff(xrange));

figure()
imagesc(real(eps));
colorbar;

neigs = 31; %modes to get
% %% mode solve
% Kx = pi/diff(xrange);
% n = 10; % num_modes;
% wvlen_min = 1;
% [Ez, Hx, Hy, omega_eigs] = solveTE_BlochX(L0, wvlen_min, xrange,...
%     yrange, eps.', Kx, Npml, n);
% 
% % for i = 1:n
% %     figure();
% %     visreal(Ez{i}, xrange, yrange);
% % end

%% MODE SCAN
band_structure = [];
wvlen_min = 1;
f = figure();
plot(linspace(0, pi/L(1)*a/(2*pi), 200), c0*1e6*linspace(0, pi/L(1)*1e-6/(2*pi*c0), 200))
hold on;
counter = 1;
for Kx = linspace(0, pi/L(1), 20)
    %Kx =0.5*pi/L(1);
    K_vec = [Kx,0];
    [Hz, Ex, Ey, omega_eigs,Am] = eigensolve_TM(L0, wvlen, xrange, ...
    yrange, eps, Npml, neigs, K_vec);
%     [Hz, Ex, Ey, omega_eigs,A] = solveTE_BlochX(L0, wvlen, xrange, yrange,...
%         eps, Kx, Npml, neigs);
%     figure(); 
%     plot(real(omega_eigs))
%     hold on;
%     plot(real(omega_eigs_mine))
    
    %% we need to filter out spurious modes (particularly modes with 
    % heavily concentrated fields near/around the PML
    filtered_eigs = [];
    for i = 1:neigs
        structure_fields = Hz{i}(:, 100:200);
        air_fields = Hz{i}(:,1:100);
        PML_fields = Hz{i}(:,1:Npml(2));
        %filter out homogeneous
        if(mean(mean(abs(PML_fields)))>1e-3)
            continue;   
        end
        if(max(max(abs(Hz{i})))-mean(mean(abs(Hz{i})))< 1e-6)
           continue; 
        end
        
        if(mean(mean(abs(structure_fields)))> mean(mean(abs(air_fields))) && ...
                real(omega_eigs(i))*1e-6/(2*pi*c0)> 0)
%             figure();
%             visreal(Hz{i}, xrange, yrange)
%             %title(num2str(i));
%             title(num2str(real(omega_eigs(i))*1e-6/(2*pi*c0)))
%             drawnow();
            filtered_eigs = [filtered_eigs, omega_eigs(i)];
        end 
    end
    disp(strcat('filtered modes: ', num2str(length(filtered_eigs))));
    band_structure = [band_structure, omega_eigs];
    for j = 1:length(omega_eigs)
        scatter(Kx*a/(2*pi),real(omega_eigs(j))*1e-6/(2*pi*c0), '.b')
        %text(Kx*a/(2*pi), real(filtered_eigs(j))*1e-6/(2*pi*c0), num2str(real(filtered_eigs(j))*1e-6/(2*pi*c0)), 'Fontsize', 10);

    end
    %plot(repmat(Kx*a/(2*pi), 1, length(filtered_eigs)), abs(filtered_eigs)*a*1e-6/(2*pi*c0), '.b');
    drawnow();
    hold on;
    counter = counter+1;
%     if(counter>1)
%         break
%     end

end

xlabel('Kx')
ylabel('eigen_frequency/(2*pi*c0)')
title('Kx band structure of the HG')
saveas(f, 'kx_band_structure_hybrid_grating.fig');



