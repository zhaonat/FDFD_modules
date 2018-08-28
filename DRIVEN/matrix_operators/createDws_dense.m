
function Dws = createDws_dense(w,s,dL,N)
    %w = 'x', 'y', or 'z'
    %s = 'b' or 'w'
    % dL: [dx dy dz] for 3D; [dx dy] for 2D
    % N: [Nx Ny Nz] for 3D; [Nx Ny] for 2D


    dw = dL('xyz' == w);  % one of dx, dy, dz
    sign = diff('bf' == s);  % +1 for s=='f'; -1 for s=='b'
    M = prod(N);  % total number of cells in domain

    ind_cur = 1:M;  % indices of current points
    ind_cur = ind_cur(:);

    ind_adj = 1:M;  % indices of adjacent (previous or next) points in the w-direction
    ind_adj = reshape(ind_adj, N);
    ind_adj = circshift(ind_adj, -sign * ('xyz' == w));
    ind_adj = ind_adj(:);

    off_diag = (sign/dw)*ones(M,1);
    on_diag  = -(sign/dw)*ones(M,1);
    Dws = sparse([ind_cur;ind_cur], [ind_adj;ind_cur], [off_diag;on_diag]);
    
%     dw = dL('xyz' == w);  % one of dx, dy, dz, all are numbers...no info about which one was selected
%     sign = diff('bf' == s);  % +1 for s=='f'; -1 for s=='b'
%     M = prod(N);  % total number of cells in domain
%     %take advantage of linear indexing
%     points = 1:M;  % indices of current points, will iterate through matrix grid.
%     points = points(:); %ensures later that points has same dims as ind_adj;
%     
%     ind_adj = 1:M;  % indices of adjacent (previous or next) points in the w-direction
%     ind_adj = reshape(ind_adj, N);%now we have a matrix of the array indices
%     ind_adj = circshift(ind_adj, -sign * ('xyz' == w)); %shift appropriately
%     ind_adj = ind_adj(:); %now what, this encodes at best about whether the offdiagonal is below or above?
%     
%     %sub2ind: matlab function to convert coordinate indices to linear index
%     linear_ind = sub2ind([M M], points, ind_adj); %points = rowsub...ind_adj = colsub
%     
%     Dws = -sign*speye(M); %M fully determines the matricial size;
%     Dws(linear_ind) = sign;
%     Dws = (1/dw)*Dws;
%     
end