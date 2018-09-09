
%% function which generates multiple unit cells on one grid

function eps = ...
    hybrid_grating_multi_unit_cell_nu(num_cells, xrange_array, yrange_array,...
    epsilon_diel,epsilon_metal, fill_factor, thickness)
    % thickness in physical units (microns)
    % fill factor (between 0 and 1)
    % we will have specified N and L for the whole simulation grid
    % beforehand
    Nx = length(xrange_array); Ny = length(yrange_array);
    Nyf = Ny/2; % these are our reference coordinates (aka origin)

    thickness_y_indices = find(abs(yrange_array) < thickness/2);
    
    eps = ones(Nx, Ny);
    eps(:, thickness_y_indices) = epsilon_diel;
    
    unit_cell_size = N(1)/num_cells; %false if x is non-uniform
    
    %% now we have to stripe it, which we will do per unit cell
    % this means, we need to locate AT LEAST, the COM of each unit cell
    % which seems to be a bad idea...
    for i = 0:num_cells-1
       x1 = i*unit_cell_size; x2 = (i+1)*unit_cell_size;
       center = (x1+x2)/2;
       xm1 = center-fill_factor*unit_cell_size/2;
       xm2 = center+fill_factor*unit_cell_size/2;
       eps(xm1:xm2,y1:y2) = epsilon_metal; 
    end
    
end

