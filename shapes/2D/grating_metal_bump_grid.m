
%% function basically is a tri-layer grating measurement
% 
function eps = grating_metal_bump_grid(N, L, epsilon_structure, fill_factor_array, ...
    thickness_array)
    %% ============================================================%%
    % epsilon_structure
    % fill_factor_array
    % thickness_array
    % L = physical dimensions of the domain
    
    %% ============================================================%%
    
    assert(length(thickness_array) == length(fill_factor_array))
    assert(sum(thickness_array) < L(2));
    Nx = N(1); Ny = N(2);

    eps = ones(N);
    num_layers = length(thickness_array);

    centerx = floor(N(1)/2);
    centery = floor(N(2)/2);
    Lx = L(1); Ly = L(2);

    % determine the overall thickness of the structure
    total_thickness = sum(thickness_array);
    %determine fractional widths of each layer
    fractional_widths = thickness_array/total_thickness;

    %determine the start boundary
    start_y = floor((Ly-total_thickness)/2*(Ny/Ly));

    %get all layer centers in y space
    for i = 1:num_layers
       y_center_array(i) = start_y+floor((thickness_array(i)/2)*(Ny/Ly));
       start_y = start_y+floor(thickness_array(i)*Ny/Ly);
    end


    for i = 1:num_layers
        x_center = centerx;
        epsilon_layer = epsilon_structure{i};
        fill_factor = fill_factor_array(i);
        epsilon_metal = epsilon_layer(1);
        epsilon_diel = epsilon_layer(2);
        thickness = thickness_array(i);

        metal_size = fill_factor*Nx;
        halfx = floor(metal_size/2);

        grating_thickness = floor(Ny*(thickness/Ly));
        halfy = floor(grating_thickness/2);

        % for each layer, we have to determine the y offset
        y_center = y_center_array(i);

        %metallic center
        eps(x_center-halfx: x_center+halfx, y_center-halfy:y_center+halfy) = epsilon_metal;
        eps(1:x_center-halfx, y_center-halfy:y_center+halfy) = epsilon_diel;
        eps(x_center+halfx+1:end, y_center-halfy:y_center+halfy) = epsilon_diel;

    end
    eps = eps.';

end

% clear
% close all
% eps = ones(100,100);
% N = size(eps);
% Nx = N(1); Ny = N(2);
% thickness_array = [0.2, 0.2, 0.2]
% L = [1,1]
% imagesc(eps)
% fill_factor_array = [0.2, 0.2, 0.2];
% epsilon_structure = {[3,6], [12,2], [3,6]};
% num_layers = length(thickness_array);
% 
% centerx = floor(N(1)/2);
% centery = floor(N(2)/2);
% Lx = L(1); Ly = L(2);
% 
% % determine the overall thickness of the structure
% total_thickness = sum(thickness_array);
% %determine fractional widths of each layer
% fractional_widths = thickness_array/total_thickness;
% 
% %determine the start boundary
% start_y = floor((Ly-total_thickness)/2*(Ny/Ly));
% 
% %get all layer centers in y space
% for i = 1:num_layers
%    y_center_array(i) = start_y+floor((thickness_array(i)/2)*(Ny/Ly));
%    start_y = start_y+floor(thickness_array(i)*Ny/Ly);
% end
% 
% 
% for i = 1:num_layers
%     x_center = centerx;
%     epsilon_layer = epsilon_structure{i};
%     fill_factor = fill_factor_array(i);
%     epsilon_metal = epsilon_layer(1);
%     epsilon_diel = epsilon_layer(2);
%     thickness = thickness_array(i);
% 
%     metal_size = fill_factor*Nx;
%     halfx = floor(metal_size/2);
% 
%     grating_thickness = floor(Ny*(thickness/Ly));
%     halfy = floor(grating_thickness/2);
% 
%     % for each layer, we have to determine the y offset
%     y_center = y_center_array(i);
% 
%     %metallic center
%     eps(x_center-halfx: x_center+halfx, y_center-halfy:y_center+halfy) = epsilon_metal;
%     eps(1:x_center-halfx, y_center-halfy:y_center+halfy) = epsilon_diel;
%     eps(x_center+halfx+1:end, y_center-halfy:y_center+halfy) = epsilon_diel;
% 
% end
% imagesc(eps)
% colorbar