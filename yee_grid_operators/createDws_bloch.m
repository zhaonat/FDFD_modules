function Dws = createDws_bloch(w, s, dL, N, k, L)
%% This function creates the sparse derivative operators in FDFD
%% with a bloch boundary condition
%% note that having a bloch boundary condition is only suitable 
%% for DRIVEN Simulation
%% for eigensolutions...you are solving a bloch eigenproblem, so the soln is different
%% Input parameters
%   w: one of 'x', 'y'
%   s: one of 'f' and 'b'
%   dL: [dx dy] for 2D
%   N: [Nx Ny] for 2D
%   k: [kx, ky] for 2D
%   L: [Lx, Ly] for 2d
%% Output parameters
%   Dws: derivative operators

%% Compute Nx, Ny
Nx = N(1); 
Ny = N(2); 

%% Sparse identity matrices
Ix = speye(Nx); 
Iy = speye(Ny); 

kx = k(1); ky = k(2);
Lx = L(1); Ly = L(2);

%% Create derivative operators
switch w
    case 'x'
        if s == 'f'
            dxf = -Ix + circshift(Ix, [0 1]);
            dxf(Nx,1) = exp(-1i*kx*Lx);
            Dws = 1/dL(1) * kron(Iy, dxf); 
        else
            dxb = Ix - circshift(Ix, [0 -1]); 
            dxb(1,Nx) = -exp(+1i*kx*Lx);
            Dws = 1/dL(1) * kron(Iy, dxb); 
        end
        
        
    case 'y'
        if s == 'f'
            dyf = -Iy + circshift(Iy, [0 1]); 
            dyf(Ny,1) = exp(-1i*ky*Ly);            
            Dws = 1/dL(2) * kron(dyf, Ix); 
        else
            dyb = Iy - circshift(Iy, [0 -1]); 
            dyb(1,Ny) = -exp(+1i*ky*Ly);            

            Dws = 1/dL(2) * kron(dyb, Ix); 
        end
end

end


%% generalized 3D version (in progress)
% the original code won't quite work because the bloch BC does not affect
% EVERY element on an off diagonal index
% moreover, depending on 'bf',
%     dw = dL('xyz' == w);  % one of dx, dy, dz
%     sign = diff('bf' == s);  % +1 for s=='f'; -1 for s=='b'
%     bloch_term = sign*exp(sign*1i*k_i*L_i);
%     M = prod(N);  % total number of cells in domain
% if(w == 'x')
%    threshold = 1;
% elseif(w =='y')
%    threshold = N(1);
% else
%    threshold = N(1)*N(2);
% end
%
% 
%     ind_cur = 1:M;  % indices of current points
%     ind_cur = ind_cur(:);
% 
%     ind_adj = 1:M;  % indices of adjacent (previous or next) points in the w-direction
%     ind_adj = reshape(ind_adj, N);
%     ind_adj_shift = circshift(ind_adj, -sign * ('xyz' == w));
%     %reflatten
%     off_diag = (sign/dw)*ones(N);
%     % the threshold depends on 'xyz' (1, Nx, Nx*Ny) respectively
%     off_diag(abs(ind_adj-ind_adj_shift) > threshold) = bloch_term;
%     off_diag = off_diag(:);
%     ind_adj = ind_adj_shift(:);
%     
%     
%     off_diag = (sign/dw)*ones(M,1);
%     % can we systematically modify off-diag to include the Bloch BC?
%     %yes, anywhere if ind_adj -circshift(ind_adj) > 1; only for x...
%
%     on_diag  = -(sign/dw)*ones(M,1);
%     Dws = sparse([ind_cur;ind_cur], [ind_adj;ind_cur], [off_diag;on_diag]);

%%


