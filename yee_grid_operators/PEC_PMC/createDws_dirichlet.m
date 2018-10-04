%% derivative operators with dirichlet boundary conditions
function Dws = createDws_dirichlet(w, s, dL,N)
    %w = 'x', 'y', or 'z'
    %s = 'b' or 'w'
    % dL: [dx dy dz] for 3D; [dx dy] for 2D
    % N: [Nx Ny Nz] for 3D; [Nx Ny] for 2D

    dw = dL('xyz' == w);  % one of dx, dy, dz, all are numbers...no info about which one was selected
    sign = diff('bf' == s);  % +1 for s=='f'; -1 for s=='b'
    M = prod(N);  % total number of cells in domain
    %take advantage of linear indexing

    ind_cur = 1:M;  % indices of current points
    ind_cur = ind_cur(:);

    ind_adj_x = 1:M;  % indices of adjacent (previous or next) points in the w-direction
    ind_adj_x = reshape(ind_adj_x, N);
    ind_adj_y = ind_adj_x;

    %% shift the row indices
    ind_adj_x = circshift(ind_adj_x, -sign * ('xyz' == w));
    ind_copy = ind_adj_x;
    %% now we need to strategically remove the indices corresponding to the PBC
    % this always corresponds to a row or column
    boundaryIndex = find(-sign*('xyz' == w));
    if(boundaryIndex == 1)
        boundaryIndex = 2;
    elseif(boundaryIndex == 2)
        boundaryIndex = 1;
    else

    end
    inds = num2cell(N);
    if(sign == -1)
        inds = repmat({1},1,length(N));
    end
    %inds = repmat({1},1,ndims(A)); 
    inds{boundaryIndex} = 1:N(boundaryIndex);
    ind_adj_x(inds{:}) = 0;
    ind_adj_y(inds{:}) = 0;
    dbc = ind_adj_x;

    ind_adj_x(ind_adj_x == 0) = [];
    ind_adj_y(ind_adj_y == 0) = [];

    xindex = [ind_adj_y.';ind_cur];
    yindex = [ind_adj_x.';ind_cur];

    off_diag = (sign/dw)*ones(length(ind_adj_x),1);
    on_diag  = -(sign/dw)*ones(M,1);
    Dws = sparse(xindex,yindex, [off_diag;on_diag]);
    Dws = Dws;

end

