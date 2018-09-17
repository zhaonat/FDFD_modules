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
            dxf(Nx,1) = exp(1i*kx*Lx);
            Dws = 1/dL(1) * kron(Iy, dxf); 
        else
            dxb = Ix - circshift(Ix, [0 -1]); 
            dxb(1,Nx) = -exp(-1i*kx*Lx);
            Dws = 1/dL(1) * kron(Iy, dxb); 
        end
        
        
    case 'y'
        if s == 'f'
            dyf = -Iy + circshift(Iy, [0 1]); 
            dyf(Ny,1) = exp(1i*ky*Ly);            
            Dws = 1/dL(2) * kron(dyf, Ix); 
        else
            dyb = Iy - circshift(Iy, [0 -1]); 
            dyb(1,Ny) = -exp(-1i*ky*Ly);            

            Dws = 1/dL(2) * kron(dyb, Ix); 
        end
end

end

