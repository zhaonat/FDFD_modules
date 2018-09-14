%% function which generates a discrete epsilon grid

function eps = add_grating_nu(eps, xrange_array, yrange_array, epsilon_edge, epsilon_center, ...
    fill_factor, thickness, y_center, x_center, lattice_constant)
    %add grating structure into an existing eps grid
    % L is typically the size of the grid, but it really specifies the unit
    % cell
    if(nargin < 8)
       x_center = L(1)/2;
    end
    if(nargin <9)
       y_center = L(2)/2; 
    end
    if(nargin <10)
       lattice_constant = L(1); 
    end
    %y_center allows us to change the height of the grating on the y axis
    % N = size of the grid in nodes
    % L has the physical dimensions

    metal_size = fill_factor*lattice_constant; %in physical units

    %% get metal bounds
    rightx = x_center+metal_size/2;
    leftx = x_center-metal_size/2;

    % get indices closest to rightx and righty
    [~,  rhalf] = min(abs(rightx-xrange_array));
    [~,  lhalf] = min(abs(leftx-xrange_array));

    %% get unit cell bounds
    lattice_bound_right =  x_center+lattice_constant/2;
    lattice_bound_left =   x_center-lattice_constant/2;

    % get indices closest to rightx and righty
    [~,  lat_rhalf] = min(abs(lattice_bound_right-xrange_array));
    [~,  lat_lhalf] = min(abs(lattice_bound_left-xrange_array));

    %% get y bounds
    upper_y = y_center +thickness/2;
    lower_y = y_center - thickness/2;
    [~,  thalf] = min(abs(upper_y-yrange_array));
    [~,  bhalf] = min(abs(lower_y-yrange_array));

    %% add grating to epsilon
    eps(lat_lhalf:lat_rhalf, bhalf:thalf) = epsilon_edge;
    eps(lat_lhalf:lat_rhalf, bhalf:thalf) = epsilon_edge;

    eps(lhalf: rhalf, bhalf:thalf) = epsilon_center;


end