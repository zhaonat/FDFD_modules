
function Dws = Laplacian9Point(w,s,dL,N)

    dw = dL('xyz' == w);  % one of dx, dy, dz, all are numbers...no info about which one was selected
    sign = diff('bf' == s);  % +1 for s=='f'; -1 for s=='b'
    M = prod(N);  % total number of cells in domain
    %take advantage of linear indexing

    ind_cur = 1:M;  % indices of current points
    ind_cur = ind_cur(:);

    ind_adj_x = 1:M;  % indices of adjacent (previous or next) points in the w-direction
    ind_adj_x = reshape(ind_adj_x, N);
    ind_adj_x0 = ind_adj_x;

    %% shift the row indices
    Dws = -sign*(5/2)*speye(M); %M fully determines the matricial size;
    for w = ['x', 'y']
        ind_adj_x = circshift(ind_adj_x0,  -sign * ('xyz' == w));
        ind_adj_2x = circshift(ind_adj_x0, -2*sign * ('xyz' == w));
        ind_adj_rx = circshift(ind_adj_x0,  sign * ('xyz' == w));
        ind_adj_r2x = circshift(ind_adj_x0, 2*sign * ('xyz' == w));

        ind_adj_xl = ind_adj_x(:);
        ind_adj_2xl = ind_adj_2x(:);
        ind_adj_rxl = ind_adj_rx(:);
        ind_adj_r2xl = ind_adj_r2x(:);

        %% conver the offdiagonals into a linear index
        linear_ind_x = sub2ind([M M], ind_cur, ind_adj_xl); %points = rowsub...ind_adj = colsub
        linear_ind_2x = sub2ind([M M], ind_cur, ind_adj_2xl); %points = rowsub...ind_adj = colsub
        linear_ind_rx = sub2ind([M M], ind_cur, ind_adj_rxl); %points = rowsub...ind_adj = colsub
        linear_ind_r2x = sub2ind([M M], ind_cur, ind_adj_r2xl); %points = rowsub...ind_adj = colsub

        Dws(linear_ind_x) = (4/3)*sign;
        Dws(linear_ind_2x) = (-1/12)*sign;
        Dws(linear_ind_rx) = (4/3)*sign;
        Dws(linear_ind_r2x) = (-1/12)*sign;

    end
    Dws = (1/dw^2)*Dws;
%     
end