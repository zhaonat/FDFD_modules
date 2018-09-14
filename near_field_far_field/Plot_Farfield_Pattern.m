function [U, Farfield] = Plot_Farfield_Pattern(r_near, norm_surf, E_near,...
    H_near, dl, theta_table, wavelength, polarization)
%function [U, Farfield] = Plot_Farfield_Patten(r_near, norm_surf, E_near, H_near, dl, theta_table, wvlen,  polarization)
%theta_table: in rad
% dl: sampling grid size
%polarization == 's' or 'p'
%U is noramlized to the average value of U. Farfield is not normalized

if polarization == 's'
    Ez_near = E_near(3,:);
    Ht_near = H_near(1:2, :);
    Ez_far_table = zeros(size(theta_table));

    for ind_theta = 1:length(theta_table)
        theta = theta_table(ind_theta);
        r_hat = [cos(theta); sin(theta)];
        Ez_far = NTFF_Transform_2D(Ez_near, Ht_near, norm_surf, r_near, dl, r_hat, wavelength, polarization);
        Ez_far_table(ind_theta) = Ez_far;
    end
    Farfield = Ez_far_table;
    
elseif polarization == 'p'
    Et_near = E_near(1:2,:);
    Hz_near = H_near(3, :);
    Hz_far_table = zeros(size(theta_table));

    for ind_theta = 1:length(theta_table)
        theta = theta_table(ind_theta);
        r_hat = [cos(theta); sin(theta)];
        Hz_far = NTFF_Transform_2D(Et_near, Hz_near, norm_surf, r_near, dl, r_hat, wavelength, polarization);
        Hz_far_table(ind_theta) = Hz_far;
    end   
    Farfield = Hz_far_table;
end
Farfield_norm = abs(Farfield)/max(abs(Farfield(:)));
U = Farfield_norm.^2;

d_theta = diff(theta_table);
d_theta = [d_theta(1), d_theta];
U = U*2*pi/sum(U(:).*d_theta(:));

%% plot
% figure;
% plot(theta_table, U, 'b', 'LineWidth', 1.5);
% xlim([0 2*pi]);
% xlabel('Theta', 'FontSize', 24);
% ylabel('Directivity', 'FontSize', 20);
% set(gca, 'FontSize', 24);

% figure;
% polar(theta_table, U, 'b');
% % xlim([0 2*pi]);
% % xlabel('Theta', 'FontSize', 24);
% % ylabel('Directivity', 'FontSize', 20);
% set(gca, 'FontSize', 24);