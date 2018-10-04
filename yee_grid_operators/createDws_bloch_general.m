function Dws = createDws_bloch_general(w, s, dL, N, k, L)

% generalized 3D version (in progress)
% the original code won't quite work because the bloch BC does not affect
% EVERY element on an off diagonal index
% moreover, depending on 'bf',
    dw = dL('xyz' == w);  % one of dx, dy, dz
    sign = diff('bf' == s);  % +1 for s=='f'; -1 for s=='b'
    bloch_term = sign*exp(sign*1i*k('xyz' == w)*L('xyz'==w));
    M = prod(N);  % total number of cells in domain
    if(w == 'x')
       threshold = 1;
    elseif(w =='y')
       threshold = N(1);
    else
       threshold = N(1)*N(2);
    end


    ind_cur = 1:M;  % indices of current points
    ind_cur = ind_cur(:);

    ind_adj = 1:M;  % indices of adjacent (previous or next) points in the w-direction
    ind_adj = reshape(ind_adj, N);
    ind_adj_shift = circshift(ind_adj, -sign * ('xyz' == w));
    %reflatten
    off_diag = (sign/dw)*ones(N);
    % the threshold depends on 'xyz' (1, Nx, Nx*Ny) respectively
    off_diag(abs(ind_adj-ind_adj_shift) > threshold) = (1/dw)*bloch_term;
    off_diag = off_diag(:);
    ind_adj = ind_adj_shift(:);
    
    % can we systematically modify off-diag to include the Bloch BC?
    %yes, anywhere if ind_adj -circshift(ind_adj) > 1; only for x...

    on_diag  = -(sign/dw)*ones(M,1);
    Dws = sparse([ind_cur;ind_cur], [ind_adj;ind_cur], [off_diag;on_diag]);

%%

end
