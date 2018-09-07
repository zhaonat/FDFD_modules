
function [Fsx, Fsy, Fsx_conj, Fsy_conj] = non_uniform_scaling(dx_scale, dy_scale)
    
    %operators which perform the row-wise scaling
    %xs: 1D array containing dx scalings (only for forward differences
    
    %% create grid of x and y points
    [Xs, Ys] = meshgrid(dx_scale, dy_scale);
    %meshgrid isn't right for y
    M = numel(Xs);

    % we have to this kind of flip because the flattening
    % operation (:) doesn't retain row-major order
    Ys=Ys'; Xs = Xs';
    Fsy = spdiags(Ys(:),0,M,M);
    Fsx = spdiags(Xs(:),0,M,M);
    
    
    % might as well construct the conjugate grid.
    xc = (dx_scale+circshift(dx_scale,[0,1]))/2;
    yc = (dy_scale+circshift(dy_scale,[0,1]))/2;
    
    [Xc, Yc] = meshgrid(xc, yc);
    Xc = Xc';
    Yc = Yc';
    Fsy_conj = spdiags(Yc(:),0,M,M);
    Fsx_conj = spdiags(Xc(:),0,M,M);
    
    
end