
function [R] = interpolate(N, w, s)  
    % this function should look VERY similar to Dws
    % w = 'x', 'y', or 'z'
    % s = 'b' or 'w'
    
    % A is a matrix with natural ordering...in fact it should be
    % diagonal...
    sign = diff('bf' == s);  % +1 for s=='f'; -1 for s=='b'
    M = prod(N);  % total number of cells in domain

    ind_cur = 1:M;  % indices of current points
    ind_cur = ind_cur(:);

    ind_adj = 1:M;  % indices of adjacent (previous or next) points in the w-direction
    ind_adj = reshape(ind_adj, N);
    ind_adj = circshift(ind_adj, -sign * ('xyz' == w));
    ind_adj = ind_adj(:);
    
    %dafuq are these points?
    off_diag = (sign/2)*ones(M,1);
    on_diag  = (sign/2)*ones(M,1);
    R = sparse([ind_cur;ind_cur], [ind_adj;ind_cur], [off_diag;on_diag]);
    
    
end