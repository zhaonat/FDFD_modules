
%% function which generates multiple unit cells on one grid

function eps_multi = ...
    hybrid_grating_multi_unit_cell_nu(num_cells, xrange_array, yrange_array,...
    epsilon_edge,epsilon_center, fill_factor, y_thickness, y_center)
    % thickness in physical units (microns)
    % fill factor (between 0 and 1)
    % we will have specified N and L for the whole simulation grid
    % beforehand
    if(nargin < 8)
       y_center = 0; 
    end
    Nx = length(xrange_array); Ny = length(yrange_array);
    eps_multi = ones(Nx, Ny);
    
    Lx = xrange_array(end)- xrange_array(1);
    lattice_constant = Lx/num_cells;
    %% now we have to stripe it, which we will do per unit cell
    % this means, we need to locate AT LEAST, the COM of each unit cell
    % which seems to be a bad idea...
    for i = 0:num_cells-1
        %determine x_center of each unit cell
        x_center = xrange_array(1)+lattice_constant*i+lattice_constant/2;
        eps_multi = add_grating_nu(eps_multi, xrange_array, yrange_array, epsilon_edge, epsilon_center, ...
    fill_factor, y_thickness, y_center, x_center, lattice_constant);

    end

    
end

