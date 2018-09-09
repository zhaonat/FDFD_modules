%% function which generates a discrete epsilon grid

function eps = hybrid_grating_grid(N, L, epsilon_diel, epsilon_metal, fill_factor, thickness)

    Nx = N(1); Ny = N(2);
    Lx = L(1); Ly = L(2);
    %initialize grid
    eps = ones(N);
    x_center = floor(Nx/2);
    y_center = floor(Ny/2);
    metal_size = fill_factor*Nx;
    halfx = floor(metal_size/2);
    
    grating_size = floor(Ny*(thickness/Ly));
    halfy = floor(grating_size/2);
    %metallic center
    eps(x_center-halfx: x_center+halfx, y_center-halfy+1:y_center+halfy) = epsilon_metal;
    eps(1:x_center-halfx, y_center-halfy+1:y_center+halfy) = epsilon_diel;
    eps(x_center+halfx+1:end, y_center-halfy+1:y_center+halfy) = epsilon_diel;


end