

function [nx,ny] = coord_to_grid(coordinate, N, xrange, yrange)
    %converts any coordinate to an equivelant index (i,j)
    % N = [Nx,Ny]
    % coordinate = [x,y]
    % xrange = [xmin, ymin]
    x = coordinate(1);
    y = coordinate(2);
    Lx = diff(xrange); Ly = diff(yrange);
    fracx = abs(x-xrange(1))/Lx;
    fracy = abs(y-yrange(1))/Ly;
    
    nx = round(fracx*N(1));
    ny = round(fracy*N(2));

    
end