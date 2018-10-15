

function [objective] = focal_objective_function(grid_field, epsilon, focus_coordinates)
    % typical objective function just requires the field distribution
    % and potentially epsilon
    % objective(E_field) = (1/2)*(E_i^2)
    % derivative of this objective w.r.t to E_i is just E_i
    % but E_i has a non-trivial derivative w.r.t to x,y and epsilon...
    % is E actually a function of the epsilon(x,y)?
    [num_coords, ~] = length(focus_coordinates);
    objective = 0;
    for i =1:num_coords
        coord = focus_coordinates(i,:);
        objective = objective + norm(grid_field(coord(1), coord(2)));
    end

end