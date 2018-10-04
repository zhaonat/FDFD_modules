


function [abc_scale] = ABC_1D(L0, Nx, dx, omega, s)
    % construct an operator for Dxf/Dxb for the Mur boundary condition
    % in 1D

    eps0 = 8.854e-12 * L0;  % vacuum permittivity in farad/L0
    mu0 = pi * 4e-7 * L0;  % vacuum permeability in henry/L0
    c0 = 1/sqrt(eps0*mu0);  % speed of light in vacuum in L0/sec
    Ix = speye(Nx);
    if s == 'f'
        dxf = -Ix + circshift(Ix, [0 1]);
        dxf(Nx,1) = 0;
        %for forward, wave is propagating right
        %dxf(Nx,Nx-1) = -1; 
        dxf(Nx,Nx) = -1+1i*omega/c0 ;
        %set the other edge to be the identity only
        dxf(1,2) =0;
        abc_scale = 1/dx * kron(1, dxf); 
    else
        dxb = Ix - circshift(Ix, [0 -1]); 
        dxb(1,Nx) = 0;
        %for the back, wave is propagating left
        dxb(1,2) = 1; dxb(1,1) = -1 - 1i*omega/c0;
        %set the other edge to be the identity only
        dxb(Nx,Nx-1) = 0;
        
        abc_scale = 1/dx * kron(1, dxb); 
        
    end
        
end