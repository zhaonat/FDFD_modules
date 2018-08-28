
function [Fsx, Fsy, Fsx_conj, Fsy_conj] = non_uniform_scaling(xs, ys)
    
    %operators which perform the row-wise scaling
    %xs: 1D array containing dx scalings (only for forward differences
    
    %% create grid of x and y points
    [Xs, Ys] = meshgrid(xs, ys);
    M = numel(Xs);
    Fsx = spdiags(Xs(:),0,M,M);
    Fsy = spdiags(Ys(:),0,M,M);
    
    
    % might as well construct the conjugate grid.
    xc = (xs+circshift(xs,-1))/2;
    yc = (ys+circshift(ys,-1))/2;
    
    [Xc, Yc] = meshgrid(xc, yc);
    Fsx_conj = spdiags(Xc(:),0,M,M);
    Fsy_conj = spdiags(Yc(:),0,M,M);
    
    
end