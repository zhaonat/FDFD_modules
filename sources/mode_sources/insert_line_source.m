

function [J] = insert_line_source(J, r_start, r_end, xrange, yrange, profile)
    % attempts to insert a particular mode profile into a J grid (Nx,Ny)
    % the source must be along a vertical or horizontal line
    % more complicated sources will be dealt separately
    % r_start = [x_start, y_start]
    % r_end = [x_end, y_end]; % point where the profile terminates;
    N = size(J);
    
    [nx0,ny0] = coord_to_grid(r_start, N, xrange, yrange);
    [nxf,nyf] = coord_to_grid(r_end, N, xrange, yrange);
    
    %% check if start and end actually forms a horizontal/vertical line
    assert (nx0-nxf == 0 | nyf -ny0 == 0, 'not a vertical or horizontal line')
   
    
    %% insert
    J(nx0:nxf, ny0:nyf) = profile;
        
    % J is modified in placee !!
    

end