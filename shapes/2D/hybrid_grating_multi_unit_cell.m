
%% function which generates multiple unit cells on one grid

function eps = ...
    hybrid_grating_multi_unit_cell(num_cells, N, L, epsilon_diel,...
    epsilon_metal, fill_factor, thickness)
    % thickness in physical units (microns)
    % fill factor (between 0 and 1)
    % we will have specified N and L for the whole simulation grid
    % beforehand
    Nx = N(1); Ny = N(2);
    dL = L./N;
    Nyf = Ny/2; % these are our reference coordinates (aka origin)
    y1 = Nyf - (thickness/2)/dL(2);
    y2 = Nyf + (thickness/2)/dL(2);
    eps = ones(N);
    eps(:, y1:y2) = epsilon_diel;
    
    unit_cell_size = N(1)/num_cells;
    
    %% now we have to stripe it, which we will do per unit cell
    for i = 0:num_cells-1
       x1 = i*unit_cell_size; x2 = (i+1)*unit_cell_size;
       center = (x1+x2)/2;
       xm1 = center-fill_factor*unit_cell_size/2;
       xm2 = center+fill_factor*unit_cell_size/2;
       eps(xm1:xm2,y1:y2) = epsilon_metal; 
    end
    
end

