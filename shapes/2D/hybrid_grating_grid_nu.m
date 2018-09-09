%% function which generates a discrete epsilon grid

function eps = hybrid_grating_grid_nu(xrange_array, yrange_array, ...
    epsilon_diel, epsilon_metal, fill_factor, thickness)

    Nx = length(xrange_array);
    Ny = length(yrange_array);
    
    %initialize grid
    eps = ones(N);
    x_center = floor(Nx/2);
    
    metal_size = fill_factor*Nx; %wrong if there is non-uniformity in x
    
    halfx = floor(metal_size/2);
    
    % grating is centered and xrange_array and yrange_array should be
    % centered
    
    grating_indices_y = find(abs(yrange_array) < thickness/2);
    
    %metallic center
    eps(x_center-halfx: x_center+halfx, grating_indices_y) = epsilon_metal;
    eps(1:x_center-halfx, grating_indices_y) = epsilon_diel;
    eps(x_center+halfx+1:end, grating_indices_y) = epsilon_diel;


end