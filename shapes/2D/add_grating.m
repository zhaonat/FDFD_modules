%% function which generates a discrete epsilon grid

function eps = add_grating(eps, L, epsilon_diel, epsilon_metal, ...
    fill_factor, thickness, y_center, x_center, lattice_constant)
    %add grating structure into an existing eps grid
    % L is typically the size of the grid, but it really specifies the unit
    % cell
    % should we account for the PML? My intuition is no...just let the
    % propagating state hit the PML and decay into it. 
    % but what if the state propagating into the PML has some imaginary
    % parts lying around.
    if(nargin < 8)
       x_center = L(1)/2;
    end
    if(nargin <7)
       y_center = L(2)/2; 
    end
    if(nargin <9)
       lattice_constant = L(1); 
    end
    %y_center allows us to change the height of the grating on the y axis
    % N = size of the grid in nodes
    % L has the physical dimensions
    N = size(eps);
    Nx = N(1); Ny = N(2);
    Lx = L(1); Ly = L(2);
    lattice_n = (lattice_constant/L(1))*N(1);

    nxcenter = round(Nx*x_center/Lx);
    nycenter = round(Ny*y_center/Ly);

    metal_size = fill_factor*lattice_n;
    halfx = floor(metal_size/2);
    lat_half = round(lattice_n/2);
    
    grating_size = floor(Ny*(thickness/Ly));
    halfy = floor(grating_size/2);
    %metallic center
    eps(nxcenter-halfx: nxcenter+halfx, nycenter-halfy:nycenter+halfy) = epsilon_metal;
    eps(nxcenter-lat_half+1:nxcenter-halfx, nycenter-halfy:nycenter+halfy) = epsilon_diel;
    eps(nxcenter+halfx+1:nxcenter+lat_half, nycenter-halfy:nycenter+halfy) = epsilon_diel;

    % sizing enforcement
    eps = eps(1:Nx, 1:Ny);

end