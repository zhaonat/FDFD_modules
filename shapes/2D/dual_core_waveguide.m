
function [eps] =  ...
    dual_core_waveguide(num_cells, N, L, epsilon_diel,...
    epsilon_metal, fill_factor, thicknesses, y_wall_1, y_wall_2, y_center_wall)
    %thicknesses: goes in order of wall1, wall2, center wall
    eps = ones(N);
    eps = hybrid_grating_multi_unit_cell_add(eps,num_cells, N, L, epsilon_diel,...
        epsilon_metal, fill_factor, thicknesses(1), y_wall_1);
    eps = hybrid_grating_multi_unit_cell_add(eps,num_cells, N, L, epsilon_diel,...
        epsilon_metal, fill_factor, thicknesses(2), y_wall_2);
    eps = hybrid_grating_multi_unit_cell_add(eps,num_cells, N, L, epsilon_diel,...
        epsilon_metal, fill_factor, thicknesses(3), y_center_wall);
end

