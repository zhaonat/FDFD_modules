%% function which generates a discrete epsilon grid

function eps = hybrid_grating_grid(N, L, epsilon_diel, epsilon_metal, ...
    fill_factor, thickness, y_center)
    
    Nx = N(1); Ny = N(2);
    Lx = L(1); Ly = L(2);
    if(nargin < 7)
        nycenter = floor(Ny/2); 
    else
        nycenter = floor(Ny*y_center/Ly);
    end
    %y_center allows us to change the height of the grating on the y axis
    % N = size of the grid in nodes
    % L has the physical dimensions

    %initialize grid
    eps = ones(N);
%     x_center = floor(Nx/2);
%     y_center = floor(Ny/2);

    nxcenter = floor(Nx/2);

    metal_size = fill_factor*Nx;
    halfx = floor(metal_size/2);
    
    grating_size = floor(Ny*(thickness/Ly));
    halfy = floor(grating_size/2);
    %metallic center
    eps(nxcenter-halfx: nxcenter+halfx, nycenter-halfy+1:nycenter+halfy) = epsilon_metal;
    eps(1:nxcenter-halfx, nycenter-halfy+1:nycenter+halfy) = epsilon_diel;
    eps(nxcenter+halfx+1:end, nycenter-halfy+1:nycenter+halfy) = epsilon_diel;


end