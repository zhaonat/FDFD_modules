


function [outputArg1,outputArg2] = epsilon_mapping(epsilon_uniform,xrange_array, ...
    yrange_array)

    %EPSILON_MAPPING Summary of this function goes here
    % We want a way to map any epsilon we define on a non-uniform grid
    % to an epsilon a uniform grid or vice versa
    % the only similar thing between uniform and non-uniform is the total 
    % domain size.
    % so what we want to get right is the location and size of the
    % structure
    % ALTERNATIVELY: given an index ix, iy on the non-uniform grid, we want
    % to know where it occurs in xrange_array and yrange_array
    
    Nu = size(epsilon_uniform);
    

end

