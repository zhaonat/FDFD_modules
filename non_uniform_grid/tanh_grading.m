
function [grading] = tanh_grading(start_value, end_value, num_points)
    
    % a tanh function appears to be smoother at the edges of the grid
    % so I'm going to try this

    center = min(start_value, end_value)+abs(end_value-start_value)/2;
    xspace = linspace(start_value, end_value, num_points);
    grading = tanh(xspace-center);
    
    
end