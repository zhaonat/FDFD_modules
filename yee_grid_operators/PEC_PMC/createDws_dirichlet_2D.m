function Dws = createDws_dirichlet_2D(w, s, dL,N)

    %% more transparent way to construct the dws matrices...but only for 2D

    Nx = N(1); Ny = N(2);
    Ix = speye(Nx); 
    Iy = speye(Ny); 

    %% Create derivative operators
    switch w
        case 'x'
            if s == 'f'
                dxf = -Ix + circshift(Ix, [0 1]);
                dxf(Nx,1) = 0;
                Dws = 1/dL(1) * kron(Iy, dxf); 
            else
                dxb = Ix - circshift(Ix, [0 -1]); 
                dxb(1,Nx) = 0;
                Dws = 1/dL(1) * kron(Iy, dxb); 
            end


        case 'y'
            if s == 'f'
                dyf = -Iy + circshift(Iy, [0 1]); 
                dyf(Ny,1) = 0;            
                Dws = 1/dL(2) * kron(dyf, Ix); 
            else
                dyb = Iy - circshift(Iy, [0 -1]); 
                dyb(1,Ny) = 0;            

                Dws = 1/dL(2) * kron(dyb, Ix); 
            end
    end   
end