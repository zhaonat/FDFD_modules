
%% function which generates multiple unit cells on one grid

function [eps,ny_bounds, metal_bound_array] = ...
    hybrid_grating_multi_unit_cell_add(eps,num_cells, N, L, epsilon_diel,...
    epsilon_metal, fill_factor, thickness, y_center)
    % thickness in physical units (microns)
    % fill factor (between 0 and 1)
    % we will have specified N and L for the whole simulation grid
    % beforehand
    lattice_constant = L(1)/num_cells;
    %metal_thickness = fill_factor*lattice_constant;
    if(nargin < 8)
        y_center = L(2)/2;
    end
    Nx = N(1); Ny = N(2);
    ny_center = round(Ny*y_center/L(2));

    dL = L./N;
    y1 = ny_center - (thickness/2)/dL(2);
    y2 = ny_center + (thickness/2)/dL(2);
    eps(:, y1:y2-1) = epsilon_diel;
    
    %unit_cell_size = N(1)/num_cells;
    metal_bound_array = [];

    %% now we have to stripe it, which we will do per unit cell
    % in fact, I wonder if it is possible to use the original grating
    % function to do it...yes with x_center..
    for i = 0:num_cells-1
%        x1 = i*unit_cell_size; x2 = (i+1)*unit_cell_size;
%        center = (x1+x2)/2;
%        xm1 = center-fill_factor*unit_cell_size/2;
%        xm2 = center+fill_factor*unit_cell_size/2;
       %eps(xm1:xm2,y1:y2) = epsilon_metal; 
       x_center = (i)*lattice_constant+(lattice_constant/2);
       [eps,ny_bounds, metalbounds] = add_grating(eps, L, epsilon_diel, epsilon_metal, ...
            fill_factor, thickness, y_center, x_center, lattice_constant);
       metal_bound_array = [metal_bound_array; metalbounds];

    end
    
end

