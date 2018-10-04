function Dws = Laplacian9pointNN(dL,N)
    
%take advantage of linear indexing
    dw = dL(1);  % one of dx, dy, dz, all are numbers...no info about which one was selected
    sign = 1;  % +1 for s=='f'; -1 for s=='b'
    M = prod(N);  % total number of cells in domain
    %take advantage of linear indexing

    ind_cur = 1:M;  % indices of current points
    ind_cur = ind_cur(:);

    ind_adj_x = 1:M;  % indices of adjacent (previous or next) points in the w-direction
    ind_adj_x = reshape(ind_adj_x, N);

    %% shift the row indices
    Dws = -sign*speye(M); %M fully determines the matricial size;
    for w = ['x', 'y']
        ind_adj_r0x = circshift(ind_adj_x, -sign * ('xy' == w));
        ind_adj_l0x = circshift(ind_adj_x, sign * ('xy' == w));
        ind_adj_rtx = circshift(ind_adj_x, sign * [1, 1]);
        ind_adj_rbx = circshift(ind_adj_x, sign * [1, -1]);
        ind_adj_ltx = circshift(ind_adj_x, sign * [-1, 1]);
        ind_adj_lbx = circshift(ind_adj_x, sign * [-1, -1]);

        ind_adj_r0x = ind_adj_r0x(:);
        ind_adj_l0x = ind_adj_l0x(:);
        ind_adj_rtx = ind_adj_rtx(:);
        ind_adj_ltx = ind_adj_ltx(:);
        ind_adj_rbx = ind_adj_rbx(:);
        ind_adj_lbx = ind_adj_lbx(:);

        %% conver the offdiagonals into a linear index
        linear_ind_r0x = sub2ind([M M], ind_cur, ind_adj_r0x); 
        linear_ind_l0x = sub2ind([M M], ind_cur, ind_adj_l0x); 
        linear_ind_rtx = sub2ind([M M], ind_cur, ind_adj_rtx); 
        linear_ind_ltx = sub2ind([M M], ind_cur, ind_adj_ltx); 
        linear_ind_rbx = sub2ind([M M], ind_cur, ind_adj_rbx); 
        linear_ind_lbx = sub2ind([M M], ind_cur, ind_adj_lbx); 

        Dws(linear_ind_r0x) = sign;
        Dws(linear_ind_l0x) = sign;
        Dws(linear_ind_rtx) = sign;
        Dws(linear_ind_ltx) = sign;
        Dws(linear_ind_rbx) = sign;
        Dws(linear_ind_lbx) = sign;

    end
    Dws = -(1/dw^2)*Dws; %% WE NEED TO DO A SIGN SWITCH FOR CORRECT ANSWER
end