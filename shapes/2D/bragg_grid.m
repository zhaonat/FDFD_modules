
% generate a bragg mirror of n layers

function [eps] = bragg_grid(N,L, H, y_offset, n_layers, epsilon, d1, d2)
    %d1: thickness of layer 1
    %d2: thickness of layer 2
    %H : [hx, hy] x dimension of the stack, y dimension of the stack
    % the stack will be centered at the bottom the grid
    

    dL = L./N;
    hx= H(1);  dx = dL(1); dy = dL(2);
    Nx = N(1);
    
    %length of the stack in terms of grid points;
    nhx = round(hx/dx);
    
    x_start = (Nx-nhx)/2+1; x_end = nhx;
    % get x_start and x_end
    ny_offset = round(y_offset/dy);
    
    % get thicknesses of each layer in the bragg stack in terms of grid
    % pointsd
    nd1 = round((d1)/dy);
    nd2 = round((d2)/dy);
    nd = [nd1,nd2];
    n_lay = nd1+nd2;
    %initialize the grid
    eps = ones(N);
    
    % start constructing the layers
    for i = 1:n_layers
        y_start =ny_offset+(i-1)*n_lay+1;
        for j = 1:2
            y_end = y_start+nd(j);
            eps(x_start:x_end, y_start:y_end) = epsilon(j);
            y_start = y_end;
        end
        
    end
    
end