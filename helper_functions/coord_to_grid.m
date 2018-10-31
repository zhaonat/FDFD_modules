

function [nx,ny] = coord_to_grid(coordinate, N, xrange, yrange)
    %converts any coordinate to an equivelant index (i,j) on the
    %discretized grid
    % param: N = [Nx,Ny]
    % param: coordinate = [x,y]
    % param: xrange = [xmin, ymax]
    % param: yrange = [ymin, ymax]

    x = coordinate(1);
    y = coordinate(2);
    Lx = diff(xrange); Ly = diff(yrange);
    fracx = abs(x-xrange(1))/Lx;
    fracy = abs(y-yrange(1))/Ly;
    
    nx = ceil(fracx*N(1));
    ny = ceil(fracy*N(2));
    if(nx == 0)
       nx = 1; 
    end
    if(ny == 0)
       ny = 1; 
    end

    
end